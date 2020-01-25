python setup.py sdist;
twine upload dist/*;
pip uninstall flyeye;
pip install flyeye;
