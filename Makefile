
test:
	py.test

sdist:
	python setup.py sdist

clean:
	rm -f *pyc; rm -f test/*pyc

realclean:	clean
	rm -rf dist
	rm -f MANIFEST
	find . -name __pycache__ | xargs rm -rf
	rm -rf .cache

