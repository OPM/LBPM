Dependencies for LBPM documentation

# install sphinx
```python 
pip install Sphinx
```
# foamatting requires sphinx read-the-docs-theme
```python
pip install sphinx-rtd-theme
```


# equation rendering requires latex and dvipng command
```
sudo apt-get install dvipng
sudo apt-get install texlive texstudio
sudo apt-get install texlive-latex-recommended texlive-pictures texlive-latex-extra
```


# To build the docs
```
Step 1) install dependencies listed above
Step 2) type 'make html' from the docs/ directory
Step 3) point your browser at ~/local/doc/build/html/index.html
```
# 
