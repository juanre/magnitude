try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup


with open("README", "r") as fp:
    long_description = fp.read()


setup(name='magnitude',
      version='0.9.4',
      description='Python library for computing with numbers with units',
      long_description=long_description,
      author='Juan Reyero',
      author_email='juan@juanreyero.com',
      url='http://juanreyero.com/open/magnitude/',
      py_modules=['magnitude'],
      classifiers=[
          "Programming Language :: Python :: 2",
          #"Programming Language :: Python :: 3",
          "Intended Audience :: Developers",
          "Intended Audience :: Science/Research",
          "License :: OSI Approved :: Apache Software License",
          "Topic :: Scientific/Engineering",
          "Topic :: Software Development :: Libraries"
      ])
