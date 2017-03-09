# Written by JP Janet for HJK Group
# Dpt of Chemical Engineering, MIT
import os, sys, copy
import glob, re, math, random, string, numpy, pybel
from math import pi
from molSimplify.Scripts.geometry import *
from molSimplify.Classes.atom3D import *
from molSimplify.Classes.mol3D import*
from molSimplify.Classes.globalvars import globalvars
from molSimplify.Informatics.graph_analyze import *
from molSimplify.Informatics.autocorrelation import *
from molSimplify.Informatics.misc_descriptors import *
def get_descriptors():
	# returns allowed descriptors
	avail_descriptors  = ['kier','truncated_kier','randic','truncated_randic',
				'con_atom_type','con_atom_en','max_con_atom_delta_en',
				'autocorrelation']
	return(avail_descriptors)
def check_descriptors(args):
	clean_descriptors = list()
	allowed_descriptors = get_descriptors()
	if args.descriptors: ## check
#		 if the descriptors are 
			     ## correctly given
		if args.debug:	
			print('descriptors are' + str(args.descriptors))
		for descriptors in args.descriptors:
			if descriptors in allowed_descriptors:
				clean_descriptors.append(descriptors)
			else:
				print('ignoring invalid descriptor ' + str(descriptors))
				print('allowed values are  '+ str(allowed_descriptors))
	else:
		clean_descriptors = allowed_descriptors
	n_desc = len(clean_descriptors)
	if args.debug:
		print('using ' + str(n_desc) + ' descriptors: ' + str(clean_descriptors) )
	if n_desc > 0:
		return clean_descriptors
	else:
		print('Error, no valid descriptors! using default instead')
		return(allowed_descriptors)
def analysis_supervisor(args):
	## main routine to control
	## analysis process
	descriptors = check_descriptors(args)
	print('more to come!')


#from sklearn import linear_model,preprocessing,metrics,feature_selection,model_selection


def get_misc_ligand_descriptors(mol):
	results_ax =  list()	
	result_eq = list()
	colnames = []
	liglist,ligdents,ligcons = ligand_breakdown(this_mol)
	ax_ligand,eq_ligand,ax_natoms,eq_natoms,ax_con,eq_con,built_ligand_list=ligand_assign(this_mol,liglist,ligdents,ligcons)
	ax_en = get_lig_EN(ax_ligand.mol,ax_con) 
	eq_en = get_lig_EN(ax_ligand.mol,eq_con)
	ax_kier = kier(ax_ligand.mol,ax_con) 
	eq_kier = kier(eq_ligand.mol,eq_con)
	ax_t_kier = kier(ax_ligand.mol,ax_con) 
	eq_t_kier = kier(eq_ligand.mol,eq_con)

	## get full ligand AC
def accquire_file(path):
	### this function reads in values from a correctly formated path
	### [name],[y],[folder_name]
	file_dict = dict()
	fail_dict = dict()
	counter = 0 # number of added paths
	ncounter = 0 # number of skipped paths 	
	if os.path.exists(path):
		with open(path,'r') as f:
			### expects csv fomart,
			### value | path
			for lines in f:
				ll = lines.split(",")
				name = ll[0]
				y_value = ll[0]
				this_path =ll[1].strip() ## check if path exists:
				
				#print(this_path)
				this_run = dft_observation(counter,this_path)
				this_run.sety(y_value)
				if os.path.isfile(this_path):
					#print('path exists')
					this_run.obtain_mol3d()
					if this_run.health:
						counter += 1
						file_dict.update({counter:this_run})
					else: ### bad geo
						this_run.comments(' geo is not healthy, culling ' +  str(this_path))
						ncounter +=1 
						fail_dict.update({counter:this_run})
				else: ### no geo file found
					this_run.comments(' geo could not be found: ' +  str(this_path))
					ncounter +=1 
					fail_dict.update({counter:this_run})
	if counter >0:
		print('file import successful, ' + str(counter) + ' geos loaded')
	return(file_dict,fail_dict)


def correlation_supervisor():
	file_dict, fail_dict = accquire_file(this_file)
 
## loop over sucessful imports to get descriptors:
#big_mat = list()

#for i,keyv in enumerate(file_dict.keys()):
	
	
	#file_dict[keyv].get_descriptor_vector(loud=False)
	#if i == 0:
		#col_names = file_dict[keyv].descriptor_names
	##print('*******************')
	#this_row = list()
	#this_row.append(float(file_dict[keyv].yvalue))
	#this_row.extend(file_dict[keyv].descriptors)
	##print(this_row)
	#big_mat.append(this_row)

#big_mat = np.array(big_mat)



##### let's do some regression
### standardize model:
#col_array = np.array(col_names)
#n_tot = len(col_array)

#X = big_mat[:,1:]
#n_obs = len(X[:,1])
#Scaler = preprocessing.StandardScaler().fit(X)
#Xs = Scaler.transform(X)

#Y = big_mat[:,0]
#Reg = linear_model.LinearRegression()
#Reg.fit(Xs,Y)
#Ypred_all_all = Reg.predict(Xs)
#rs_all_all = metrics.r2_score(Y,Ypred_all_all)
##print('the full R2 is  '+str(rs_all))

#loo = model_selection.LeaveOneOut()
#r_reduce = list()
#mse_reduce = list()
### stepwise reduce the feature set until only one is left
#for n in range(0,n_tot):
	#reductor = feature_selection.RFE(Reg,n_tot-n,step=1,verbose=0)
	#reductor.fit(Xs,Y)
	#Ypred_all = reductor.predict(Xs)
	#rs_all = metrics.r2_score(Y,Ypred_all)
	#mse_all = metrics.mean_squared_error(Y,Ypred_all)
	#r_reduce.append(rs_all)
	#mse_reduce.append(mse_all)
### reduce to one feature

#reductor_features = list()
#for i,ranks in enumerate(reductor.ranking_):
	#reductor_features.append([col_array[i], ranks])
#reductor_features =  sorted(reductor_features,key=lambda x: x[1] )
##print(reductor_features)

#print('****************************************')

### select best number using cv
#selector = feature_selection.RFECV(Reg,step=1,cv=loo,verbose=0,scoring='neg_mean_squared_error')
#selector.fit(Xs,Y)
#select_mse = selector.grid_scores_ 
#Ypred = selector.predict(Xs)

#rs = metrics.r2_score(Y,Ypred)
#n_opt = selector.n_features_
#print('analzyed ' +  str(n_obs) +  ' molecules')
#print('optimal number of features is ' + str(n_opt) + ' of total ' + str(n_tot))
#print('the full-space R2 is  '+str(rs_all_all))
#print('the opt R2 is  '+str(rs))

#opt_features = col_array[selector.support_]

#ranked_features = list()
#for i,ranks in enumerate(selector.ranking_):
	#ranked_features.append([col_array[i], ranks])
#ranked_features =  sorted(ranked_features,key=lambda x: x[1] )
##print(ranked_features)


#X_r = selector.transform(Xs)

#reg_red = linear_model.LinearRegression()
#reg_red.fit(X_r,Y)
#Ypred_r = reg_red.predict(X_r)
#errors = [ Y[i] - Ypred_r[i] for i in range(0,n_obs) ]
#coefs = reg_red.coef_
#intercept = reg_red.intercept_
#mse_all = metrics.mean_squared_error(Y,Ypred_all_all)
#mse_r = metrics.mean_squared_error(Y,Ypred_r)
#print('the optimal variables are: ' + str(opt_features))
#print('the coefficients are' + str(coefs))
#print('the intercept is '+ str(intercept))
#print('tthe MSE is '+ str(mse_r))
#print('the MSE (ALL) is '+ str(mse_all))

### train
#with open('rankings.csv','w') as f:
	#f.write('RFE_rank,RFE_col,RFECV_rank,RFECV_col, \n')
	#for i,items in enumerate(reductor_features):
		#f.write(str(items[0])  + ',' + str(items[1])+ ',' +
			#str(ranked_features[i][0])  + ','+ str(ranked_features[i][1])+ '\n')



#with open('y_data.csv','w') as f:
	#for items in Y:
		#f.write(str(items) + '\n')
#with open('y_pred_r.csv','w') as f:
	#for items in Ypred_r:
		#f.write(str(items) + '\n')
#with open('y_full_all.csv','w') as f:
	#for items in Ypred_all_all:
		#f.write(str(items) + '\n')
#with open('rfe_r.csv','w') as f:
	#for items in r_reduce:
		#f.write(str(items) + '\n')
#with open('rfe_mse.csv','w') as f:
	#for items in mse_reduce:
		#f.write(str(items) + '\n')
#with open('select_mse.csv','w') as f:
	#for items in select_mse:
		#f.write(str(items) + '\n')
#with open('errors.csv','w') as f:
	#for items in errors:
		#f.write(str(items) + '\n')
#with open('transformed_space.csv','w') as f:
	#for i in range(0,n_obs):
		#for j in range(0,n_opt):
			#if j == (n_opt-1):
				#f.write(str(X_r[i][j])+'\n')
			#else:
				#f.write(str(X_r[i][j])+',')
#with open('full_space.csv','w') as f:
	#for names in col_names:
		#f.write(names+',')
	#f.write('\n')
	#for i in range(0,n_obs):
		#for j in range(0,n_tot):
			#if j == (n_tot-1):
				#f.write(str(Xs[i][j])+'\n')
			#else:
				#f.write(str(Xs[i][j])+',')
#with open('scaling.csv','w') as f:
	#means = Scaler.mean_
	#var = Scaler.var_
	#f.write('name, mean,variance \n')
	#for i in range(0,n_tot):
		#f.write(str(col_names[i])+','+str(means[i]) +','+str(var[i])+','+str(selector.ranking_[i])+'\n')
