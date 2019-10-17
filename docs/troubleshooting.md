# Troubleshooting Bactopia
It was bound to happen, an error has occurred or a bug has shown itself. Now let's see if we can fix it!

!!! error "Don't see your error/bug? Post an issue on GitHub"
    If you've encountered an error or bug not seen here, please post an issue at [Bactopia GitHub Issues](https://github.com/bactopia/bactopia/issues). This will help greatly to track down the error and fix it!

## Past Errors

!!! error "These should have been fixed."
    The following errors should have been fixed. If you are still receiving them, please make sure you are using the most up to date version of Bactopia. If you are, and you are still recieving one of these errors, please reopen the associated issue on GitHub. Thanks!

### Failed to create Conda Environment

!!! info "This was fixed in v1.2.1"
    A new function `bactopia build` was introduced that builds the Conda environments outside of Nextflow. If you are still receiving this error please reopen the [Better handling of conda environments?](https://github.com/bactopia/bactopia/issues/13) issue.

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
