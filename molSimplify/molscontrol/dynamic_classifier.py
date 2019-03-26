import os
import time
import logging
from pkg_resources import resource_filename, Requirement
from collections import OrderedDict
import numpy as np
import pickle
import json
from keras.models import load_model
from io_tools import obtain_jobinfo, read_geometry_to_mol, get_geo_metrics, get_bond_order, get_gradient, \
    get_mullcharge, kill_job
from clf_tools import get_layer_outputs, dist_neighbor, get_entropy


class dft_control:
    '''
    Attribute of the class:
    step_now: current step for an optimization
    self.mode: mode of the dynamic classifier: "full": geo metrics + electronic structure as descriptors, stable for
    terachem users, and "geo": only geometry metrics as descriptors, which can be used for all kinds of quantum
    chesmitry software.
    self.mode_allowed: allowed modes.
    self.step_decisions: steps at which the dynamic classifier can make predictions. (Resizing to be implemented).
    self.scrpath: path to the scratch directory.
    self.geofile: filename of the optimization trajectory in a xyz file format. This is the minimum requirement to use 
    the dynamic classifier.
    self.bofile: filename of the trajectory for the bond order matrix (for mode = "full").
    self.chargefile: filename of the trajectory for the Mulliken charge(for mode = "full").
    self.gradfile: filename of the trajectory for the gradient matrix (for mode = "full").
    self.mols: a list of mol3D objects.
    self.features: a combinaton of all feaures of {"feature_name1": [var_step0, var_step1, ...]}
    self.features_norm: self.features after normalization.
    self.feature_mat: input for the dynamic classifier. Dim: (1, step_now, number_features).
    self.preditions: a dictionary of predictions (probabilities of success) with the format of {step_now: prediction}
    self.lses: similar as above but for the latent space entropy (LSE).
    self.train_data: Training (data, labels) for the dynamic classifiers.
    self.status: status of a simulation. True for live and false for dead.
    self.modelfile: path to the model file.
    self.models: a dictory of dynamic classifiers with the format of {step_now: model}
    self.normalization_dict: a dictionary for data normalizations. For geometry metrics, we use value/cutoff
    as standarization, and for other descriptors, we use (value-mean_step0)/std_step0 as standardization.
    self.features_dict: a dictionary for features.
    self.avrg_latent_dist_train: averaged 5-neighbor distance in the latent space for the training data.
    self.lse_cutoff: cutoff for the LSE of model confidence on the prediction.
    self.debug: Whether in debug mode. True means testing on a complete set of files with a finished job. False means
    on-the-fly job control.
    self.pid: pid to kill.
    '''

    def __init__(self, mode='full',
                 scrpath='./scr/',
                 geofile='optim.xyz',
                 bofile='bond_order.list',
                 chargefile='charge_mull.xls',
                 gradfile='grad.xyz',
                 modelsfile=False,
                 normfile=False,
                 traindatafile=False,
                 dataname=False,
                 initxyzfile="initgeo.xyz",
                 modelname='conv',
                 normvecname='norm_dict.json',
                 logfile='./molscontrol.log',
                 lse_cutoff=0.3,
                 debug=False,
                 pid=False):
        self.step_now = -1
        self.mode = mode
        self.mode_allowed = ["full", "geo"]
        self.step_decisions = [2, 5, 10, 15, 20, 30, 40]
        self.scrpath = scrpath
        self.initxyzfile = initxyzfile
        self.geofile = geofile
        self.bofile = bofile
        self.chargefile = chargefile
        self.gradfile = gradfile
        self.mols = list()
        self.features = OrderedDict()
        self.features_norm = OrderedDict()
        self.feature_mat = list()
        self.preditions = OrderedDict()
        self.lses = OrderedDict()
        self.train_data = list()
        self.status = True
        self.modelfile = modelsfile
        self.normfile = normfile
        self.traindatafile = traindatafile
        self.modelname = modelname
        self.dataname = dataname
        self.normvecname = normvecname
        self.logfile = logfile
        self.models = dict()
        self.normalization_dict = dict()
        self.lse_cutoff = lse_cutoff
        self.debug = debug
        self.pid = pid
        self.features_dict = {"full": {0: 'bo_0', 1: 'bo_sv0', 2: 'bo_offsv0', 3: 'bo_sv1', 4: 'bo_offsv1', 5: 'bo_sv2',
                                       6: 'bo_offsv2', 7: 'bo_sv3', 8: 'bo_offsv3', 9: 'bo_eq_mean', 10: 'bo_ax_mean',
                                       11: 'grad_0', 12: 'grad_sv0', 13: 'grad_intsv0', 14: 'grad_sv1',
                                       15: 'grad_intsv1',
                                       16: 'grad_sv2', 17: 'grad_intsv2', 18: 'grad_maxnorm', 19: 'grad_intmaxnorm',
                                       20: 'grad_rms', 21: 'grad_eq_mean', 22: 'grad_ax_mean', 23: 'charge_0',
                                       24: 'charge_eq_mean', 25: 'charge_ax_mean', 26: 'flag_oct',
                                       27: 'inspect_oct_angle_devi_max', 28: 'inspect_max_del_sig_angle',
                                       29: 'inspect_dist_del_all', 30: 'inspect_dist_del_eq',
                                       31: 'inspect_devi_linear_avrg',
                                       32: 'inspect_devi_linear_max', 33: 'actural_rmsd_max'},
                              "geo": {0: 'flag_oct', 1: 'inspect_oct_angle_devi_max', 2: 'inspect_max_del_sig_angle',
                                      3: 'inspect_dist_del_all', 4: 'inspect_dist_del_eq',
                                      5: 'inspect_devi_linear_avrg',
                                      6: 'inspect_devi_linear_max', 7: 'actural_rmsd_max'}
                              }
        self.avrg_latent_dist_train = {"full": {2: 6.34, 5: 7.59, 10: 4.83, 15: 5.21,
                                                20: 5.06, 30: 9.34, 40: 8.70},
                                       "geo": {2: 4.0876, 5: 5.3512, 10: 5.8573, 15: 6.1089,
                                               20: 6.4067, 30: 10.4459, 40: 9.4269}
                                       }
        self.files_track = {"full": {self.geofile: 0, self.bofile: 0, self.gradfile: 0, self.chargefile: 0},
                            "geo": {self.geofile: 0}
                            }
        self.file_updated = {"full": {self.geofile: False, self.bofile: False,
                                      self.gradfile: False, self.chargefile: False},
                             "geo": {self.geofile: False}
                             }
        self.init_mol = read_geometry_to_mol(self.initxyzfile)
        self.job_info = obtain_jobinfo(self.initxyzfile)
        self.initialize_logger()
        self.initialize_features()
        self.load_models()
        self.load_training_data()
        self.load_normalization_vec()
        self.initilize_file_track_dict()

    def initialize_logger(self):
        logging.basicConfig(filename=self.logfile, filemode='w', level=logging.DEBUG,
                            format='%(name)s - %(levelname)s - %(message)s')
        logging.info('----Logger for the dynamic classifier for on-the-fly job control---')
        logging.info('Monitoring job with PID %s' % str(self.pid))
        if not self.pid:
            logging.warning('NO PID is inputed. Cannot do cany control.')

    def initialize_features(self):
        try:
            for idx, fname in self.features_dict[self.mode].items():
                self.features.update({fname: []})
            logging.info('Feature initialized.')
        except Exception as e:
            logging.error('Feature initialization failed.', exc_info=True)

    def get_file_path(self, filein):
        return self.scrpath + '/' + filein

    def load_models(self):
        # modelpath = self.modelfile + '/' + self.mode
        if not self.modelfile:
            modelpath = resource_filename(Requirement.parse("molSimplify"),
                                          "molSimplify/molscontrol/models/" + self.mode + "/")
        else:
            modelpath = self.modelfile
            logging.warning("Using user-specified models from %s." % modelpath)
        try:
            for step in self.step_decisions:
                _model = '/%s_%d.h5' % (self.modelname, step)
                _modelname = modelpath + _model
                logging.info("Loading model: %s ..." % _modelname.split('/')[-1])
                self.models.update({step: load_model(_modelname)})
        except Exception as e:
            logging.error('Failed at model loading.', exc_info=True)

    def load_training_data(self):
        # datapath = self.traindatafile + '/' + self.mode + '/' + self.dataname
        if not self.traindatafile:
            datapath = resource_filename(Requirement.parse("molSimplify"),
                                         "molSimplify/molscontrol/data/" + self.mode + "/train_data.pkl")
        else:
            datapath = self.traindatafile
            logging.warning("Using user-specified models from %s." % datapath)
        try:
            with open(datapath, 'rb') as f:
                _train_data = pickle.load(f)
            for key, val in _train_data.items():
                self.train_data.append(val)
            logging.info("Training data loaded.")
        except Exception as e:
            logging.error('Failed at training data loading.', exc_info=True)

    def load_normalization_vec(self):
        # normvecpath = self.normfile + '/' + self.mode + '/' + self.normvecname
        if not self.normfile:
            normvecpath = resource_filename(Requirement.parse("molSimplify"),
                                            "molSimplify/molscontrol/normalization_vec/" + self.mode + "/norm_dict.json")
        else:
            normvecpath = self.normfile
            logging.warning("Using user-specified models from %s." % normvecpath)
        try:
            with open(normvecpath, 'rb') as f:
                self.normalization_dict = json.load(f)
            logging.info('Normalization vectors loaded')
        except:
            logging.error('Failed at normalization vector loading.', exc_info=True)

    def update_features(self):
        dict_combined = {}
        if self.mode in self.mode_allowed:
            geometrics_dict = get_geo_metrics(init_mol=self.init_mol, job_info=self.job_info,
                                              geofile=self.get_file_path(self.geofile))
            dict_combined.update(geometrics_dict)
        if self.mode == 'full':
            bondorder_dict = get_bond_order(bofile=self.get_file_path(self.bofile),
                                            job_info=self.job_info, num_sv=4)
            dict_combined.update(bondorder_dict)
            gradient_dict = get_gradient(gradfile=self.get_file_path(self.gradfile),
                                         job_info=self.job_info, num_sv=3)
            dict_combined.update(gradient_dict)
            mullcharge_dict = get_mullcharge(chargefile=self.get_file_path(self.chargefile),
                                             job_info=self.job_info)
            dict_combined.update(mullcharge_dict)
        for idx, fname in self.features_dict[self.mode].items():
            self.features[fname].append(dict_combined[fname])

    def normalize_features(self):
        for idx, fname in self.features_dict[self.mode].items():
            if fname in self.features_dict["geo"].values():
                self.features_norm[fname] = np.array(self.features[fname]) / self.normalization_dict[fname]
            else:
                self.features_norm[fname] = (np.array(self.features[fname]) - self.normalization_dict[fname][0]) / \
                                            self.normalization_dict[fname][1]

    def prepare_feature_mat(self):
        self.feature_mat = list()
        for fname, vec in self.features_norm.items():
            self.feature_mat.append(vec)
        self.feature_mat = np.transpose(self.feature_mat)
        self.feature_mat = self.feature_mat.reshape(1, self.step_now + 1, len(self.features.keys()))

    def predict(self):
        pred = self.models[self.step_now].predict(self.feature_mat)
        self.preditions.update({self.step_now: pred[0][0]})
        logging.info('Prediction made at step %d.' % self.step_now)
        logging.info('Prediction summary, %s' % str(self.preditions))

    def calculate_lse(self):
        train_latent = get_layer_outputs(self.models[self.step_now], -4, self.train_data[0][:, :self.step_now + 1, :],
                                         training_flag=False)
        test_latent = get_layer_outputs(self.models[self.step_now], -4, self.feature_mat, training_flag=False)
        __, nn_dists, nn_labels = dist_neighbor(test_latent, train_latent, self.train_data[1],
                                                l=5, dist_ref=self.avrg_latent_dist_train[self.mode][self.step_now])
        lse = get_entropy(nn_dists, nn_labels)
        self.lses.update({self.step_now: lse[0]})
        logging.info('LSE summary, %s' % str(self.lses))

    def make_decision(self):
        if self.preditions[self.step_now] <= 0.5 and self.lses[self.step_now] < self.lse_cutoff:
            logging.critical("!!!!job killed at step %d!!!!!!" % self.step_now)
            logging.info("Reasons: a prediction of %.4f with LSE of %.4f" % (self.preditions[self.step_now],
                                                                             self.lses[self.step_now]))
            self.status = False
            kill_job(self.pid)
        else:
            logging.info('This calculation seems good for now at step %d' % self.step_now)

    def initilize_file_track_dict(self):
        existed = False
        logging.info("Initialize the files to be tracked during the geometry optimization.")
        logging.info("This may take a while until the first step of SCF calculation finishes...")
        while not existed:
            for filename in self.files_track[self.mode]:
                filepath = self.get_file_path(filename)
                if os.path.isfile(filepath):
                    self.file_updated[self.mode].update({filename: True})
            existed = all(value == True for value in self.file_updated[self.mode].values())
        for filename in self.files_track[self.mode]:
            filepath = self.get_file_path(filename)
            self.files_track[self.mode].update({filename: os.path.getmtime(filepath)})
            self.file_updated[self.mode].update({filename: False})
        logging.info("Tracking files initialization completes.")
        time.sleep(3)
        ### Gather features for the 0th step of the optimization
        try:
            self.update_features()
            self.normalize_features()
            self.step_now += 1
            self.prepare_feature_mat()
            logging.info("%d step feature obatined." % self.step_now)
        except Exception as e:
            logging.warning('Cannot obtain the information of the zeroth step.', exc_info=True)

    def check_updates(self):
        for filename, val in self.files_track[self.mode].items():
            filepath = self.get_file_path(filename)
            if os.path.getmtime(filepath) - val > 3:
                self.file_updated[self.mode].update({filename: True})
        updated = all(value == True for value in self.file_updated[self.mode].values())
        if self.debug:
            updated = True
        if updated:
            for filename in self.files_track[self.mode]:
                filepath = self.get_file_path(filename)
                self.files_track[self.mode].update({filename: os.path.getmtime(filepath)})
                self.file_updated[self.mode].update({filename: False})
            self.step_now += 1
        return updated

    def update_and_predict(self):
        stop = False
        updated = self.check_updates()
        if self.debug:
            updated = True
        if self.status:
            if updated:
                self.update_features()
                self.normalize_features()
                self.prepare_feature_mat()
                if self.step_now in self.step_decisions:
                    self.predict()
                    self.calculate_lse()
                    self.make_decision()
                else:
                    logging.info("At step %d, decision is not activated." % self.step_now)
            if self.step_now > max(self.step_decisions):
                logging.warning(
                    "Step number is larger than the maximum step number that we can make dicision (%d steps). The dynamic classifer is now in principle deactivated." % max(
                        self.step_decisions))
                stop = True
        else:
            stop = True
        return stop
