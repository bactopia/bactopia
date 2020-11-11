# Update Docs
Currently have a janky set up for docs. The Markdown files for the docs are on the `bactopia` repo but the docs are built to `bactopia.github.io` repo.

To setup local `mkdocs` envrionment:
```
screen -S bactopia-docs
conda create -n mkdocs -c conda-forge 
conda activate mkdocs 
pip install mkdocs-material
pip install mkdocs-material-extensions
pip install pygments
pip install pymdown-extensions
git pull git@github.com:bactopia/bactopia.git
git pull git@github.com:bactopia/bactopia.github.io.git
```

In the `bactopia` repo directory, run: `mkdocs serve`

This will serve it on `http://0.0.0.0:8000`.

Once everything is good and you are ready to push an update to the docs you will need to change directories:

```
cd ../bactopia.github.io/
mkdocs gh-deploy -f ../bactopia/mkdocs.yml -b master
```

This will push the docs to the `master` branch on `bactopia.github.io`.

To serve the docs again locally:
```
cd ../bactopia
mkdocs serve
```

Now you can exit the screen (CTRL+A CTRL+D)
