
python -m pydoc -w photovoltaic > help.html
python setup.py sdist
python setup.py bdist_wheel

pip uninstall photovoltaic
pip install dist\photovoltaic-0.1.8-py3-none-any.whl
pause