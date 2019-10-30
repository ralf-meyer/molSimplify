from setuptools import setup, find_packages

setup(name='molSimplify', version="v1.2.7-alpha", packages=find_packages(),
      entry_points={'console_scripts': ['molsimplify = molSimplify.__main__:main',
                                        'molscontrol = molSimplify.molscontrol.molscontrol:main']},
      package_data={
          'molsimplify': ["Data/*.dat", "Data/*.dat", "Bind/*.dat", "Ligands/*.dict", "icons/*.png"]
      },
      data_files=[("molSimplify", ["molSimplify/Data/ML.dat"])],
      install_requirements=['numpy','scipy','scikit-learn','pandas','keras==2.2.5','tensorflow==1.4.0'],
      setup_requires=['pytest-runner'],
      tests_require=['pytest'],
      include_package_data=True
      )
