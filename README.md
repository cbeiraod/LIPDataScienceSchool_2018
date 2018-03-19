# LIPDataScienceSchool_2018
Data challenge from the LIP School "Data Science in (Astro)particle physics and the Bridge to Industry"


## Create virtualenv and install stuff

Create it:
```sh
virtualenv --system-site-packages -p [pythonVer] keras-tf
```

Install stuff:
```sh
source keras-tf/bin/activate
pip install --upgrade tensorflow keras pandas scikit-learn root_numpy h5py matplotlib
pip install --upgrade git+git://github.com/albertbup/deep-belief-network.git
```

Once finished:
```sh
deactivate
```
