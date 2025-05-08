from setuptools import setup, find_namespace_packages


setup(
    name='plasmid-surgery',
    version='0.1.0',
    packages=find_namespace_packages(),
    entry_points={
        "console_scripts": [
            "plasmid-surgery = plasmid_surgery.__main__:main",
        ]
    },
    scripts=[],
    package_data={
    },
    install_requires=[
    ],
    description='',
    url='https://github.com/BCCDC-PHL/plasmid-surgery',
    author='Dan Fornika',
    author_email='dan.fornika@bccdc.ca',
    include_package_data=True,
    keywords=[],
    zip_safe=False
)
