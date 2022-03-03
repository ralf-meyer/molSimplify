from setuptools import setup, find_packages

setup(name='molSimplify',
      version='v1.6.0',
      packages=find_packages(),
      entry_points={
          'console_scripts': [
              'molsimplify = molSimplify.__main__:main',
              'molscontrol = molSimplify.molscontrol.molscontrol:main',
<<<<<<< HEAD
              'jobmanager = molSimplify.job_manager.resub:main']
      },
      package_dir={'molSimplify': 'molSimplify'},
      package_data={
          'molSimplify': ['Data/*.dat', 'Bind/*.dat', 'Ligands/*.dict',
                          'icons/*.png', 'python_nn/*.csv', 'python_krr/*.csv',
                          'tf_nn/*/*', 'molscontrol/*/*']
      },
      data_files=[('molSimplify', ['molSimplify/Data/ML.dat'])],
      install_requires=['numpy', 'scipy', 'scikit-learn',
                        'pandas', 'keras', 'tensorflow', 'pyyaml'],
      setup_requires=['pytest-runner'],  # this may result some package conflict in local conda build. comment it out if needed.
=======
              'jobmanager = molSimplify.job_manager.resub:main',
              'molsimplify-optimize = molSimplify.optimize.main:main']
      },
      package_dir={'molSimplify': 'molSimplify'},
      package_data={
          'molSimplify': ["Data/*.dat", "Bind/*.dat", "Ligands/*.dict",
                          "icons/*.png", "python_nn/*.csv", "python_krr/*.csv",
                          "tf_nn/*/*", "molscontrol/*/*"]
      },
      data_files=[("molSimplify", ["molSimplify/Data/ML.dat"])],
      install_requires=['numpy', 'scipy', 'scikit-learn', 'pandas', 'keras',
                        'tensorflow', 'pyyaml', 'ase'],
      setup_requires=['pytest-runner'], # this may result some package conflict in local conda build. comment it out if needed.
<<<<<<< HEAD
>>>>>>> e11021d (Add molsimplify-optimize cmdl script)
      tests_require=['pytest'],
=======
      tests_require=['pytest', 'numdifftools'],
>>>>>>> 674c6d0 (Custom implementation of numerical_hessian)
      include_package_data=True
      )
