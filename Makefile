
sdist:
	python setup.py sdist

clean:
	rm -f *pyc; rm -f test/*pyc

realclean:	clean
	rm -rf dist
	rm -f MANIFEST

