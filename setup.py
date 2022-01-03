from setuptools import setup, find_packages

with open("README.md", mode="r", encoding="utf-8") as readme_file:
    readme = readme_file.read()



setup(
      name="SCProcessing",
      version="0.0.2",
      author="A. Ali Heydari",
      author_email="aliheydari@ucdavis.edu",
      description="Preprocessing single cell RNA-seq data for machine learning purposes ",
      long_description=readme,
      long_description_content_type="text/markdown",
      license="MIT",
      url="https://github.com/dr-aheydari/SCRealVAE",
      download_url="https://github.com/SindiLab/SCData-PreProcessing",
      packages=find_packages(),
      install_requires=[
                        'scanpy',
                        'anndata',
                        'pandas',
                        'scipy',
                        'numpy',
                        'tqdm',
                        'python-igraph',
                        #'fa2' only needed for older version of scanpy (< 0.6.1),
                        'leidenalg'
                        ],
      classifiers=[
                   "Development Status :: 2 - Beta",
                   "Intended Audience :: Science/Research",
                   "License :: OSI Approved :: MIT Software License",
                   "Programming Language :: Python :: 3.6",
                   "Topic :: Scientific/Engineering :: Genomics :: Bioinformatics :: Machine Learning"
                   ],
      keywords="Single Cell RNA-seq, Data Preprocessing, Sparse Single Cell, Annotated Data"
      )
