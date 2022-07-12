from setuptools import setup,find_packages
import os


__VERSION__ = '0.2.12'
#BASE_DIR = os.path.dirname(os.path.abspath(__file__))  ##
BASE_DIR = os.path.dirname(__file__)
os.chdir(BASE_DIR)

ext_modules = []
try:
    from pybind11.setup_helpers import Pybind11Extension 
    import pybind11
    ext_modules=[Pybind11Extension(name='PyAGH.FCOEFF',  # 模块名称
                            sources=['scrc/pybind11_fcoeff.cpp'],    # 源码
                            define_macros = [('VERSION_INFO', __VERSION__)],
                            language='c++',
                            cxx_std=11,
                            include_dirs = [pybind11.get_include()],
                           )]


except ImportError:
    pass


setup(

    name='PyAGH',  # package name
    version= __VERSION__,  # package version
    description='A package for calculating kinship matrix',  # package description
    author='Zhao Wei',
    author_email='852322127@qq.com',
    packages=find_packages(),
    zip_safe=False,
    python_requires='>=3.5',
    setup_requires = ["pybind11"],
    install_requires=["numpy",
                      "pandas >=1.1.0",
                      "sympy",
                      "scipy",
                      "polars",
                      "pybind11",
                      "matplotlib",
                      ],
    include_package_data=True,
    license='MIT License',
    platforms=["all"],
    classifiers=[
          'Development Status :: 3 - Alpha',
          'Natural Language :: English',
          'License :: OSI Approved :: MIT License',
          'Intended Audience :: Science/Research',
          'Operating System :: OS Independent',
          'Programming Language :: Python',
          'Programming Language :: Python :: 3',
          'Programming Language :: Python :: 3.5',
          'Programming Language :: Python :: 3.6',
          'Programming Language :: Python :: 3.7',
          'Programming Language :: Python :: 3.8',
          'Programming Language :: C++',
          'Topic :: Scientific/Engineering :: Bio-Informatics',
      ],
    ext_modules=ext_modules,
    package_data={
        '':['data/*'],
               },
    
    
)