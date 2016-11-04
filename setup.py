from setuptools import setup,find_packages

setup(name='molSimplify',version="v1.1-alpha",packages=find_packages(),
      entry_points={'console_scripts': ['molsimplify = molSimplify.__main__:main']},
      package_data={
          'molsimplify':["Data/*.dat","Data/*.dat","Bind/*.dat","Ligands/*.dict"]
      },
      include_package_data = True,
     )
