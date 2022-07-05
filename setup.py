from setuptools import setup,find_packages



VERSION = '0.1.1'

setup(

    name='PyAGH',  # package name
    version=VERSION,  # package version
    description='A package for calculating kinship matrix',  # package description
    author='Zhao Wei',
    author_email='852322127@qq.com',
    packages=find_packages(),
    zip_safe=False,
    python_requires='>=3.5',
    install_requires=["numpy",
                      "pandas",
                      "sympy",
                      "scipy",
                      "polars "],
    package_data={
        'PyAGH': ['*.so']
    },
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
    
)