python setup.py sdist;
python setup.py test;
twine upload dist/*;
pip uninstall flyeye;
pip install flyeye;
