from setuptools import find_packages, setup

setup(
    name="swipe_modules",
    version="0.1.0",
    description="Routines that adapt litebird_sim to LSPE-SWIPE",
    author="Luca Pagano",
    author_email="luca.pagano@unife.it",
    zip_safe=False,
    packages=find_packages(),
    python_requires=">=3.7.1,<3.10",
    install_requires=[
        "numpy>=1.18",
        "numba>=0.54",
        "astropy>=4.0",
        "litebird_sim",
        "setuptools",
        "mpi4py",
    ],
    include_package_data=True,
)
