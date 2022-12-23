import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="clustertools",
    version="1.0.8",
    author="Jeremy J. Webb",
    author_email="webb@astro.utoronto.ca",
    description="A python packaged for analysing star clusters",
    long_description=long_description,
    long_description_content_type='text/markdown',
    license='MIT',
    packages=["clustertools","clustertools/analysis","clustertools/cluster","clustertools/io","clustertools/io/data","clustertools/tidaltail","clustertools/util"],
    setup_requires=['numpy>1.8','scipy'],
    install_requires=['galpy','numba'],
    include_package_data=True,
    )
