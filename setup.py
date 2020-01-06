from setuptools import find_packages, setup

with open('LICENSE', encoding='utf-8') as f:
    LICENSE = f.read()

setup(
    name='logiq',
    version='0.1.3',
    description='Logiq, a simple to use quantum simulator library',
    long_description='Github repo: https://github.com/Bnz-0/logiq',
    author='Bnz',
    author_email='matteo.benzi97@gmail.com',
    url='https://github.com/Bnz-0/logiq',
    license=LICENSE,
    python_requires='>=3.6.0',
    packages=find_packages(exclude=('tests', 'docs', 'TODO')),
    install_requires=['numpy']
)
