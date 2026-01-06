from setuptools import setup, find_packages

setup(
    name="DynaPIN",
    version="0.1.1",
    packages=find_packages(include=["DynaPIN", "DynaPIN.*"]),
    python_requires=">=3.10",
    install_requires=[
        "interfacea",
    ],
    entry_points={
        "console_scripts": [
            "dynapin=DynaPIN.dynapin:main",
        ],
    },
    include_package_data=True,
)
