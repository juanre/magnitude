
test:
	py.test

sdist:
	python setup.py sdist

install:
	sudo python setup.py install

clean:
	rm -f *pyc; rm -f test/*pyc

realclean:	clean
	sudo rm -rf dist
	sudo rm -rf *.egg-info
	rm -f MANIFEST
	find . -name __pycache__ | xargs rm -rf
	rm -rf .cache

