#!/usr/bin/env python
import io
import re
from setuptools import setup, find_packages
import sys

with io.open('./DNPy/__init__.py', encoding='utf8') as version_file:
    version_match = re.search(r"^__version__ = ['\"]([^'\"]*)['\"]", version_file.read(), re.M)
    if version_match:
        version = version_match.group(1)
    else:
        raise RuntimeError("Unable to find version string.")


with io.open('README.rst', encoding='utf8') as readme:
    long_description = readme.read()


setup(
    name='DNPy',
    version=version,
    description='ODNP evaluation with an UI',
    long_description=long_description,
    author='M. Hadi Timachi',
    author_email='h.timachi@gmail.com',
    license='GNU General Public License v2 (GPLv2)',
    packages=find_packages(
        exclude=[
            'docs', 'tests',
            'windows', 'macOS', 'linux',
            'iOS', 'android',
            'django'
        ]
    ),
    classifiers=[
        'Development Status :: 1 - Planning',
        'License :: OSI Approved :: GNU General Public License v2 (GPLv2)',
    ],
    install_requires=[
    ],
    options={
        'app': {
            'formal_name': 'DNPy UI',
            'bundle': 'com.spintoolbox'
        },

        # Desktop/laptop deployments
        'macos': {
            'app_requires': [
                'toga-cocoa==0.3.0.dev8',
            ]
        },
        'linux': {
            'app_requires': [
                'toga-gtk==0.3.0.dev8',
            ]
        },
        'windows': {
            'app_requires': [
                'toga-winforms==0.3.0.dev8',
            ]
        },

        # Mobile deployments
        'ios': {
            'app_requires': [
                'toga-ios==0.3.0.dev8',
            ]
        },
        'android': {
            'app_requires': [
                'toga-android==0.3.0.dev8',
            ]
        },

        # Web deployments
        'django': {
            'app_requires': [
                'toga-django==0.3.0.dev8',
            ]
        },
    }
)
