# linsync cluster scripts

# Set up

On your cluster, you will need a few tweaks before the scripts will work for you:

1. Create sub-folders `props`, `processshellscripts`, `log`, `err`
1. `runProcess.sh` - at lines 51-60 edit the qsub command to be precisely what your cluster needs. This will include editing your memory/resource requirements if they are supplied in the qsub command to your cluster.
1. `startProcessTemplate.sh` - at lines 18-21 edit the resources you require if your cluster checks resources from the script that will be run. Also edit which Matlab or other modules you need to load in.
1. Run `chmod u+x *.sh` so that all of your shell scripts here are executable

# How it works

In brief:
1. You set up the parameters for your experiments in `parameters.m`.
1. You call `runManyProcesses.sh 1 1 200` for example to set up 200 parallel experiments -- these run several repeat sample networks each in parallel over the same parameter settings, rather than parallelising over parameters. Each job requires 3-4 hours on my cluster.
1. `runManyProcesses.sh` creates a lot of dummy copies of the parameters file (basically only changing the results filename for each) and then calls `runProcess.sh` to submit a cluster job for each. There are probably simpler ways to do this with array jobs on the cluster, but I had this running already ...
1. Each cluster job is to call `startProcessTemplate.sh`, which in turn runs Matlab with `runComputeSyncResults.m` with the given parameter file.
1. When the cluster jobs are finished, you have the 200 (or however many you ran) results files in the sub-folder of results.
1. You can call `cleanup.sh` to remove the temporary files.
1. The `combineSyncResults.m` script then pulls all of these together into a single results file; you can run this on your laptop after pulling all the results files across.

# Example - recreating Figure 1 from the paper

The parameters file `figure1ParametersCluster.m` is provided for use on the cluster, specifying only 10 repeats per job.
It is almost the same as that supplied in `../2023-AnalyticRelationshipPaper`.
Also see how a parameter argument will be substituted as `[@P1]` in the results file name, this ensures the parallel jobs all write to different results files.
So we copy that here as `parameters.m` for the cluster scripts to use:
```shell
cp ./figure1ParametersCluster.m ./parameters.m
```

Next, make sure that the directory shown in parameters.folder exists, aside from the `[@P1]` subfolder (else writing the results will fail). E.g.:
```shell
mkdir results
mkdir results/N100-randRing-d4-b1.00-c0.50-dir-k4-cont
```

Now we will submit 200 cluster jobs running `runComputeSyncResults.m` via `startProcessTemplate.sh`:
```shell
./runManyProcesses.sh 1 1 200
```

Then we wait -- each job should take 2-3 hours, but may take longer to reach the head of the queue.

Once all the jobs are finished, you should have mat files for each value of `S` in folders `1` through `200`
under your results folder.
Now we combine these results into single `.mat` files for each `S`.
You can do this on the cluster itself if you are allowed user-level Matlab jobs there, or else copy all of the 
results folders back to your own machine.
In either case, now make sure that you set `parameters.combineResultsFrom` to where all of the results folders are.
Then you can run in Matlab:
```matlab
combineSyncResults(parameters);
```
or passing the same parameters file in that you used on the cluster.

We now have single results files for each value of `S` (in `parameters.folder`) and are now ready to
plot the results in Matlab.
Before we plot them, we just need to adjust `parameters.repeats` from the number of repeats
in each thread (10 above) to the total across all threads (2000 above):

```matlab
% Adjust the parameters so that the results files are found properly:
parameters.repeats = 2000;
% Now plot the results:
plotErrorInEmpiricalSyncResults(parameters);
figure(2); % Select the correct plot we're keeping
print('-depsc', [parameters.folder, 'fig1.eps'])
saveas(gca, [parameters.folder, 'fig1.fig'], 'fig');
```

# Example - recreating Figure 4 from the paper

The parameters file `figure4aParametersCluster.m` is provided for use on the cluster, specifying only 10 repeats per job.
It is almost the same as that supplied in `../2023-AnalyticRelationshipPaper`.
Also see how a parameter argument will be substituted as `[@P1]` in the results file name, this ensures the parallel jobs all write to different results files.
So we copy that here as `parameters.m` for the cluster scripts to use:
```shell
cp ./figure4aParametersCluster.m ./parameters.m
```

Next, make sure that the directory shown in parameters.folder exists, aside from the `[@P1]` subfolder (else writing the results will fail). E.g. for 4a (change for the others):
```shell
mkdir results
mkdir results/N100-randRing-d2-b1.00-c0.50-dir-k4-cont
```

Now we will submit 200 cluster jobs running `runComputeSyncResults.m` via `startProcessTemplate.sh`:
```shell
./runManyProcesses.sh 1 1 200
```

Then we wait -- each job should take up to a few minutes, but may take longer to reach the head of the queue.

Once all the jobs are finished, you should have mat files in folders `1` through `200`
under your results folder.
Now we combine these results into a single `.mat` file.
You can do this on the cluster itself if you are allowed user-level Matlab jobs there, or else copy all of the 
results folders back to your own machine.
In either case, now make sure that you set `parameters.combineResultsFrom` to where all of the results folders are.
Then you can run in Matlab:
```matlab
combineSyncResults(parameters);
```
or passing the same parameters file in that you used on the cluster.

We now have a single result file for this parameter set (in `parameters.folder`) and are now ready to
plot the results in Matlab.
Before we plot them, we just need to adjust `parameters.repeats` from the number of repeats
in each thread (10 above) to the total across all threads (2000 above):

```matlab
% Adjust the parameters so that the results files are found properly:
parameters.repeats = 2000;
% Now plot the results:
plotSyncResults(parameters);
figure(1); % Select the correct plot we're keeping
print('-depsc', [parameters.folder, 'fig4a.eps'])
saveas(gca, [parameters.folder, 'fig4a.fig'], 'fig');
```

To reproduce subfigures 4b-f also, I've supplied sample parameters files here.
You can basically re-reun this section replacing 4a with the subfigure label for the next experiment.
