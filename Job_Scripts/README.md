# Deleterious_GP/Job_Scripts
This directory contains job scripts that run on the 
[University of Minnesota Supercomputing Institute](https://www.msi.umn.edu/)
compute clusters. The
queueing system at MSI uses PBS TORQUE.

## Task Arrays
Scripts in these directories that perform batch processing of samples use a
PBS TORQUE feature called _task arrays_. They allow users to process lists of
samples without submitting jobs with a loop. A task array and multiple job
submission with a loop perform effectively the same thing, but do not require
the user to maintain a separate script for each job.

To use a task array, use the `-t [RANGE]` option with `qsub`. Within the script,
the variable `${PBS_ARRAYID}` will take values as specified in `[RANGE]`. The
submitted script must contain some special variables for task arrays to work. 
For instance, a script called `my_script.job` may contain the following lines:

```bash
SAMPLES=($(cat /path/to/sample/list.txt))
CURRENT_SAMPLE=${SAMPLES[${PBS_ARRAYID}]}
samtools view view -f 3 -F 256 -bS ${CURRENT_SAMPLE} > ${CURRENT_SAMPLE/sam/bam}
```

The first line creates an array from a text file that contains a list of paths,
one per line (much like a [Consed](http://www.phrap.org/consed/consed.html)
file-of-files), then uses the `${PBS_ARRAYID}` variable to get a specific
file from it, and uses a samtools trimming command on it.

By running

    $ qsub -t 0-10 my_script.job

you will submit the script to the queue, and run 11 total jobs, one for each
value of `${PBS_ARRAYID}` from 0 to 10, inclusive. If you run `qstat`, you will
see `[]` after your job ID, which tells you it is an array job.

To delete jobs from a task array, you can specify the array index as part of the
`qdel` command. `qdel 12345[0]` will delete the first job out of the array. You
can supposedly specify ranges here, so `qdel 12345[0-10]` should work, too, but
I haven't seen that work on MSI.
