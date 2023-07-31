# SETUP

# import the hcmv data set
hcmv
# store the number of palindromes
PALINDROMES = 296
# store the number of bases along the DNA sequence
BASES = 229354


# PROBLEM 1: simulate 296 palindrome sites chosen at random along a DNA sequence
# of 229,354 bases using RNG

# collection of all base locations on the DNA sequence
locations = 1:BASES
# simulate a uniform distribution of 296 palindrome sites chosen at random
# 5 times
for (i in 1:5) {
  # sample 296 locations from the possible locations
  simulation = sample(locations, PALINDROMES)
  # graph a histogram of each simulation
  hist(simulation, 
       main = paste("Random Simulation #", i, sep = ''), 
       xlab = "Location")
}


# PROBLEM 2: use graphical methods to examine the spacings between consecutive
# palindromes and sum of consecutive pairs, triplets, etc, spacings

# function to calculate the consecutive sums of pairs of an array
# parameters: array as input to calculate consecutive sums of pairs for
# returns: new_array with consecutive sums of pairs of input array
consecutive_pairs = function(array) {
  new_array = array(dim = length(array) - 1)
  for (i in 1:length(new_array)) {
    new_array[i] = array[i] + array[i + 1]
  }
  return(new_array)
}

# function to calculate the consecutive sums of triplets of an array
# parameters: array as input to calculate consecutive sums of triplets for
# returns: new_array with consecutive sums of triplets of input array
consecutive_triplets = function(array) {
  new_array = array(dim = length(array) - 2)
  for (i in 1:length(new_array)) {
    new_array[i] = array[i] + array[i + 1] + array[i + 2]
  }
  return(new_array)
}

# function to calculate the consecutive differences of an array
# parameters: array as input to calculate consecutive differences for
# returns: new_array with consecutive differences of input array
consecutive_difference = function(array) {
  new_array = array(dim = length(array) - 1)
  for (i in 1:length(new_array)) {
    new_array[i] = abs(array[i] - array[i + 1])
  }
  return(new_array)
}

# convert the data to an array
data_array = hcmv[["location"]]
# find the spacings between consecutive pairs in the data
spacing_pairs = consecutive_difference(consecutive_pairs(data_array))
# find the spacings between consecutive triplets in the data
spacing_triplets = consecutive_difference(consecutive_triplets(data_array))

# compute the spacing for placing palindromes at equal distances
null_spacing = as.integer(BASES/PALINDROMES)
# create an array to store palindromes placed at equal distances
null_data = array(dim = PALINDROMES)
# fill the array with locations equally spaced from each other
for (i in 1:length(null_data)) {
  null_data[i] = i*null_spacing
}

# find the spacing of consecutive pairs for the null data
null_spacing_pairs = consecutive_difference(consecutive_pairs(null_data))
# find the spacing of consecutive triplets for the null data
null_spacing_triplets = consecutive_difference(consecutive_triplets(null_data))

# store arguments for graphing the spacing histograms as variables
spacing_xlab = "Spacing"
spacing_xlim = c(0, 200000)
spacing_ylim = c(0, 90)

# run 5 simulations
for (i in 1:5) {
  # simulation of uniform scatter for what we expect to see as an array
  random_data = sample(locations, PALINDROMES)
  # find the spacing of consecutive pairs for the null data
  random_spacing_pairs = 
    consecutive_difference(consecutive_pairs(random_data))
  # find the spacing of consecutive triplets for the null data
  random_spacing_triplets = 
    consecutive_difference(consecutive_triplets(random_data))
  # plot the simulation results
  hist(random_spacing_pairs,
       main = 
         paste("Spacing Between Consecutive Pairs (Random Simulation #", i, ')', 
               sep = ''),
       xlab = spacing_xlab,
       xlim = spacing_xlim,
       ylim = spacing_ylim)
  hist(random_spacing_triplets,
       main = 
         paste("Spacing Between Consecutive Triplets (Random Simulation #", i, 
               ')',
                sep = ''),
       xlab = spacing_xlab,
       xlim = spacing_xlim,
       ylim = spacing_ylim)
}

# graph the pair / triplet spacings as histograms
hist(spacing_pairs, 
     main = "Spacing Between Consecutive Pairs", 
     xlab = spacing_xlab,
     xlim = spacing_xlim,
     ylim = spacing_ylim)
hist(spacing_triplets,
     main = "Spacing Between Consecutive Triplets",
     xlab = spacing_xlab,
     xlim = spacing_xlim,
     ylim = spacing_ylim)
hist(null_spacing_pairs,
     main = "Spacing Between Consecutive Pairs Under the Null Distribution",
     xlab = spacing_xlab,
     xlim = spacing_xlim,
     ylim = spacing_ylim)
hist(null_spacing_triplets,
     main = "Spacing Between Consecutive Triplets Under the Null Distribution",
     xlab = spacing_xlab,
     xlim = spacing_xlim,
     ylim = spacing_ylim)


# PROBLEM 3: examine the counts of palindromes in various regions of the DNA

# function to fill in 0's for missing intervals in interval counting 
# data
# parameters: data frame to fill in 0's for, intervals integer value to fill
# in until the data frame is complete
# returns: completed data frame  
fill_counts = function (data, intervals) {
  data_intervals = data[["interval"]]
  missing_intervals = intervals - length(data_intervals)
  fill_intervals = array(dim = missing_intervals)
  index = 1
  for (i in 1:intervals) {
    if (is.element(i, data_intervals) == FALSE) {
      fill_intervals[index] = i
      index = index + 1
    }
  }
  fill_zeros = array(0, dim = missing_intervals)
  fill_data = data.frame(fill_intervals, fill_zeros)
  names(fill_data)[1] = "interval"
  names(fill_data)[2] = "count"
  data = rbind(data, fill_data)
  data = data[order(data["interval"]),]
  row.names(data) = NULL
  return(data)
}

# function to categorize data into any given number of equal length
# intervals
# parameters: data as a data frame with palindrome locations to categorize into
# intervals, intervals integer value to split the data by
# returns: data frame with column corresponding to which interval the location
# falls within
categorize_intervals = function (data, intervals) {
  new_data = data.frame(data)
  for (i in 1:intervals) {
    if (i == 1) {
      new_data["interval"] = 
        replace(new_data["location"], 
                (new_data["location"] >= 1) & 
                  (new_data["location"] <= as.integer(BASES*(1/intervals))), 1)
    }
    else {
      new_data["interval"] = 
        replace(new_data["interval"], 
                (new_data["location"] > as.integer(BASES*(i - 1)/intervals)) & 
                  (new_data["location"] <= as.integer(BASES*(i/intervals))), i)
    }
  }
  return(new_data)
}

# function to find frequency of palindrome sites in a data set for any given
# number of equally spaced intervals
# parameters: data as data frame to count frequency for, intervals integer
# value to categorize the data into
# returns: new data frame of intervals in the data with number of palindromes
# for each interval
intervals_frequency = function (data, intervals) {
  data = categorize_intervals(data, intervals)
  data_counts = aggregate(data["interval"], by = data["interval"], FUN = length)
  names(data_counts)[2] = "count"
  data_counts = fill_counts(data_counts, intervals)
  return(data_counts)
}

# declare the number of intervals to equally split the data by
intervals = c(10, 25, 50, 75, 100)
# sort the intervals array by least to greatest intervals value
intervals = intervals[order(intervals)]
# count how many interval values to run simulations for
intervals_count = length(intervals)
# create an empty array to store observed TVDs
tvds_observed = array(NA, dim = intervals_count)
# create an empty array to store TVD p-values
tvd_p_values = array(NA, dim = intervals_count)
# create an empty array to store Chi-Squared p-values
chisq_p_values = array(NA, dim = intervals_count)
# create an empty array to store the highest frequency of any interval
max_counts = array(NA, dim = intervals_count)
# create an empty array to store the highest frequency of any interval as
# a proportion
max_proportions = array(NA, dim = intervals_count)
# create an empty array to store the interval with the greatest number of
# palindrome sites
max_intervals = array(NA, dim = intervals_count)

# iterate through each interval value
for (x in 1:intervals_count) {
  counts = intervals_frequency(hcmv, intervals[x])
  # find the highest number of locations in any one interval
  max_count = max(counts["count"])
  # store the highest frequency in the array of highest frequencies
  max_counts[x] = max_count
  # find the interval with the greatest number of palindromes
  max_interval = subset(counts, counts["count"] == max_count)[["interval"]]
  # store the interval with the highest frequency of palindromes in the array
  # of intervals with highest frequency of palindromes
  max_intervals[x] = max_interval
  # find the proportions of the data that fall into each interval
  proportions = counts["count"]/PALINDROMES
  # find the highest proportion of the locations in any one interval
  max_proportion = max(proportions["count"])
  # store the highest proportion in the array of highest proportions
  max_proportions[x] = max_proportion
  # create the theoretical null distribution
  null = array(data = 1/intervals[x], dim = intervals[x])
  # calculate the observed TVD for the data
  tvd_observed = sum(abs(proportions - null))
  # store the observed TVD in the array for observed TVDs
  tvds_observed[x] = tvd_observed
  
  # declare the number of simulations to run under the null hypothesis
  simulations = 500
  
  # create an array to store the simulation results in
  results = array(NA, dim = simulations)
  # run the simulations
  for (j in 1:simulations) {
    # simulate a distribution under the null hypothesis
    simulation = sample(locations, PALINDROMES)
    # create a data frame using the randomized data
    simulation = data.frame(simulation)
    # rename the data frame's column to "location"
    names(simulation)[1] = "location"
    simulation_counts = intervals_frequency(simulation, intervals[x])
    # find the proportions of the randomized data that fall into each interval
    simulation_proportions = simulation_counts["count"]/PALINDROMES
    # calculate the simulated TVD for the trial
    tvd_simulated = sum(abs(simulation_proportions - null))
    # store the TVD statistic in the results array
    results[j] = tvd_simulated
  }
  # calculate the p-value
  p_value = mean(results >= tvd_observed)
  # store the p-value in the array of TVD p-values
  tvd_p_values[x] = p_value
  # plot the random distribution of the TVD
  hist(results, 
       main = 
         paste("Random Distribution of the TVD (", intervals[x], " Intervals)", 
               sep = ''), 
       xlab = "Simulation TVDs")
  # add the observed value as a point on the histogram
  points(tvd_observed, 0, col = "red", pch = 19)
  
  # perform Pearson's Chi-Squared Test on the frequency of locations for each
  # interval
  chisq_results = chisq.test(counts[["count"]])
  # fetch the p-value for the Chi-Squared Test
  chisq_p_value = chisq_results$p.value
  # store the p-value in the list of Chi-Squared p-values
  chisq_p_values[x] = chisq_p_value
}
# compile the simulation results into a data frame
simulation_results = 
  data.frame(intervals, max_intervals, max_counts, max_proportions, 
             tvds_observed, tvd_p_values, chisq_p_values)


# PROBLEM 4: analyze the interval with the greatest number of palindromes

# create an empty array to store the start of each biggest cluster
start = array(NA, dim = intervals_count)
# create an empty array to store the end of each biggest cluster
end = array(NA, dim = intervals_count)
# iterate through the intervals
for (i in 1:intervals_count) {
  # fetch the row corresponding to the current interval value
  interval_entry = 
    subset(simulation_results, 
           simulation_results["intervals"] == intervals[i])
  # find the cluster with the highest frequency of palindromes
  cluster = interval_entry[["max_intervals"]]
  # if the cluster is the first cluster in the data
  if (cluster == 1) {
    # then the start of this cluster is location 1
    start[i] = 1
  }
  # otherwise this cluster is not the first cluster in the data
  else {
    # so calculate the start of this cluster
    start[i] = as.integer(BASES*((cluster - 1)/intervals[i]))
  }
  # compute the end of this cluster
  end[i] = as.integer(BASES*(cluster/intervals[i]))
}
# compile the interval starts and ends into a data frame
interval_data = data.frame(intervals, start, end)
# merge the start and end columns into the simulation results
simulation_results = merge(simulation_results, interval_data)
# iterate through the intervals
for (i in 1:intervals_count) {
  # first interval to compare
  interval1 = intervals[i]
  interval1_entry = 
    subset(simulation_results, 
           simulation_results["intervals"] == interval1)
  # first interval's start
  interval1_start = interval1_entry[["start"]]
  # first interval's end
  interval1_end = interval1_entry[["end"]]
  # create an empty array to store if interval2 is in interval1
  inside = array(NA, dim = intervals_count)
  # iterate through the intervals
  for (j in 1:intervals_count) {
    # second interval to compare
    interval2 = intervals[j]
    # if the intervals are the same, or if the first interval is wider than
    # the second interval
    if ((i == j) | (interval1 > interval2)) {
      # encode this entry as NA
      inside[j] = NA
    }
    # otherwise, the intervals are different, and the first interval is shorter
    # than the second interval
    else {
      # so fetch the entry for the second interval
      interval2_entry = 
        subset(simulation_results, 
               simulation_results["intervals"] == interval2)
      # second interval's start
      interval2_start = interval2_entry[["start"]]
      # second interval's end
      interval2_end = interval2_entry[["end"]]
      # if the second interval is contained within the first interval
      if ((interval2_start >= interval1_start) & 
          (interval2_end <= interval1_end)) {
        # encode this entry as TRUE
        inside[j] = TRUE
      }
      # otherwise, the second interval is not inside the first interval
      else {
        # so encode this entry as FALSE
        inside[j] = FALSE
      }
    }
  }
  # create a data frame using the recorded boolean values
  inside_data = data.frame(intervals, inside)
  # rename the new column as "inside_X"
  names(inside_data)[2] = paste("inside_", interval1, sep = '')
  # join the simulation results with this data
  simulation_results = merge(simulation_results, inside_data, all.x = TRUE)
}


# ADVANCED ANALYSIS: quantify the difference between spacings distributions

# declare the number of simulations to run
simulations = 500
# create an empty array to store KS statistics in for pair spacings
pairs_results = array(NA, dim = simulations)
# create an empty array to store KS statistics in for triplet spacings
triplets_results = array(NA, dim = simulations)
# run the simulations
for (i in 1:simulations) {
  # randomly sample 296 locations from the possible locations
  simulation = sample(locations, PALINDROMES)
  # find the spacings between consecutive sum of pairs
  simulation_spacing_pairs = 
    consecutive_difference(consecutive_pairs(simulation))
  # find the spacings between consecutive sum of triplets
  simulation_spacing_triplets = 
    consecutive_difference(consecutive_triplets(simulation))
  # store the KS test results for pair spacings
  ks_result_pairs = ks.test(spacing_pairs, simulation_spacing_pairs)
  # store the KS test results for triplet spacings
  ks_result_triplets = ks.test(spacing_triplets, simulation_spacing_triplets)
  # obtain the KS test statistic for pair spacings
  ks_result_pairs_diff = ks_result_pairs$statistic
  # obtain the KS test statistic for triplet spacings
  ks_result_triplets_diff = ks_result_triplets$statistic
  # store the KS test statistic for pair spacings
  pairs_results[i] = ks_result_pairs_diff
  # store the KS test statistic for triplet spacings
  triplets_results[i] = ks_result_triplets_diff
}
# find the average KS test statistic for pair spacings
mean_ks_pairs = mean(pairs_results)
# find the average KS test statistic for triplet spacings
mean_ks_triplets = mean(triplets_results)
