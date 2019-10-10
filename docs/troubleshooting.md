# Troubleshooting Bactopia
It was bound to happen, an error has occurred or a bug has shown itself. Now let's see if we can fix it!

!!! error "Don't see your error/bug? Post an issue on GitHub"
    If you've encountered an error or bug not seen here, please post an issue at [Bactopia GitHub Issues](https://github.com/bactopia/bactopia/issues). This will help greatly to track down the error and fix it!

### Failed to create Conda Environment
Occasionally on the first run of Bactopia you will encounter this error: 
```
Caused by:
  Failed to create Conda environment
  command: conda env create --prefix /data/apps/bactopia/conda/cache/envs/bactopia-gather_fastqs-8bddd22dc63ec39a71c4ea8fd7704f7a --file /data/apps/bactopia/conda/gather_fastqs.yml
  status : 1
  message:
    InvalidArchiveError("Error with archive /home/rpetit/miniconda3/envs/bactopia/pkgs/python-3.7.3-h33d41f4_1.tar.bz2.  You probably need to delete and re-download or re-create this file.  Message from libarchive was:\n\nFailed to create dir 'lib'",)
```

While it may look like this is related to Nextflow, it is actually a Conda error that occurs when installing multiple environments at once which Nextflow likes to do. This can also occur from a timeout or loss of internet connectivity (under a different error name).

**Recommended Solution**  
Using `--dry_run` at runtime will run Bactopia with dummy data on a single core to prevent parallel creation of conda environments. This process will take only as long as it takes to create the environments.

If you have suggestions for how to better handle this, check out the submitted [Better handling of conda environments?](https://github.com/bactopia/bactopia/issues/13) issue. Your feedback would be greatly appreciated!

