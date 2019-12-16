#SUPPLEMENTARY MATERIAL for
#Modelling capacity along a patient pathway with delays to transfer and discharge
#Authors: Richard Wood, Ben Murch
#Date: June 2018

################## #
#1. USER INPUTS ####
################## #

wd<-"H:/My Documents/"      #location where "input_file.csv" is located

nruns<- 10000              #the number of runs to perform within the simulation
sim_time<- 365              #the number of time units (e.g. days) to simulate
front_door_queue <- TRUE    #choose whether to allow a queue to form at the front door (TRUE) or reject arrivals if at capacity (FALSE)

############### #
############### #

#CONTENTS

#1. USER INPUTS (line 7)
#2. SYSTEM SETUP (line 57)
#3. READ IN DATA FROM CSV (line 77)
#4. THE SIMULATION FUNCTION (line 104)
  #4.1 Destination-choosing function
  #4.2 Simulation results matrix setup
    #4.2.1 Column names
    #4.2.2 Initial values
  #4.3 Event parameters
  #4.4 Simulation run
    #4.4.1 Event list
      #4.4.1.1 arrivals
      #4.4.1.2 active treatment/service completion
      #4.4.1.3 delayed discharge completion
    #4.4.2 Action upon event
      #4.4.2.1 arrivals
      #4.4.2.2 active treatment/service completion
        #4.4.2.2.1 single phase
        #4.4.2.2.2 multi-phase
      #4.4.2.3 Delayed discharge completion
    #4.4.3 Upstream transfers
    #4.4.4 Occupancy counts
    #4.4.5 Checks
  #4.5 Outputs (single simulation run)
    #4.5.1 Simulation matrix
    #4.5.2 Plots
      #4.5.2.1 Plot 1 (occupancy and blocked)
      #4.5.2.2 Plot 2 (occupancy and queues)
    #4.5.3 Performance results (summary statistics)
      #4.5.3.1 Time at occupancy/blockage level calculations
      #4.5.3.2 Summary statistics matrix
#5. PARALLELISATION AND ITERATION (line 1065)
#6. AGGREGATION OF RESULTS ACROSS SIMULATION RUNS (line 1088)



################### #
#2. SYSTEM SETUP ####
################### #

setwd(wd)

wu_time<- 100        #the number of time units (e.g. days) to use as a warm-up period
init_occ<- 0.85       #the initial occupancy of the unit at the start of the simulation
seedstart<- 542       #the random number seed to use

##install.packages("doSNOW")
require(plyr)
require(foreach)
require(doSNOW)

#create new folder for output results
folder_name<-paste0("/output_", nruns, "runs_", format(Sys.time(), "%Y%m%d_%H%M%S"))
dir.create(paste0(getwd(),folder_name))


############################ #
#3. READ IN DATA FROM CSV ####
############################ #

dat<-read.csv(paste0(getwd(),"/input_file.csv"),row.names=1)

#tr is the transition matrix - it's the data from the input file without the columns for bed numbers, arrivals or LOS
tr<-dat[,1:nrow(dat)]
units<-rownames(dat)
trans<-rownames(dat)[which(!is.na(dat$beds))]
final<-rownames(dat)[which(is.na(dat$beds))]
finaldd<-rownames(dat)[which(!is.na(dat$dd))]

ar<-list()
beds<-list()
los<-list()
dd<-list()
num_phases <- list()
dd_phases <- list()

for (i in trans) ar[[i]]<-dat$ar[which(rownames(dat)==i)]
for (i in trans) beds[[i]]<-dat$beds[which(rownames(dat)==i)]
for (i in trans) los[[i]]<-dat$los[which(rownames(dat)==i)]
for (i in finaldd) dd[[i]]<-dat$dd[which(rownames(dat)==i)]
for (i in trans) num_phases[[i]] <- dat$phases[which(rownames(dat) == i)]
for (i in finaldd) dd_phases[[i]] <- dat$ddphases[which(rownames(dat) == i)]


############################## #
#4. THE SIMULATION FUNCTION ####
############################## #

sim_run<-function() {
  
  ##################################### #
  #4.1 Destination-choosing function ####
  
  #This function is called in section B.2., "Active Service Completion", of the code below.
  #It is used to choose which transfer/discharge destination to send a patient who has completed the last phase of their active stay in a
  #transient unit to (where the input "x" is named vector of cumulative sums of the probabilities of discharge to the given units/destinations
  #and the input "y" is a single sample from the uniform distribution)
  fn1<-function(x,y) {
    for (i in 1:length(x)) if (y<x[i]) break
    return(names(x)[i])
  }
  
  ######################################## #
  #4.2 Siumulation results matrix setup ####
  
  #automated setup of simulation matrix "res"
  ##first create a vector "n" of column names for the matrix
  
  ###################### #
  #4.2.1 Column names ####
  
  # "t" - event time (total time elapsed up to this point - distinct from "time" below, which records how long the system has remained in a given state)
  # it records cumulative event time in decimal days
  
  n<-"t"
  for (i in trans) {
    #creates columns for the "front door" queue for units which receive arrivals from outside the system
    if (is.na(ar[[i]])==FALSE) n<-c(n,paste0(i,"_q"))
    #creates columns to count the total number (sum of patients in each phase) of patients currently receiving active treatment at a unit
    n<-c(n,paste0(i,"_act"))
    #creates occupancy columns for each active LOS phase, for units which have multi-phase LOS (one column per phase)
    #(single-phase units/services have no phase columns, just the "act" column)
    if (num_phases[[i]] > 1)
      for (k in 1:num_phases[[i]]) {
        n<- c(n, paste0(i, "_phase-", k))
      }
    #identifies all destinations to which patients can be transferred or discharged
    dch_dest<-colnames(tr)[which(!is.na(tr[which(rownames(tr)==i),]))]
    #identifies transient states (units which are transferred to within the system)
    dch_dest1<-unlist(intersect(dch_dest,trans))
    #creates columns for internal queues to count delayed transfers
    #(a separate queue for each pair of units/services patients can be transferred between prior to discharge)
    if (length(dch_dest1)>0) n<-c(n,paste0(i,"_dc_",dch_dest1))
    #identifies discharge destinations (absorbing states - outside the system) which have a delay parameter
    final_delays <- unlist(intersect(dch_dest,finaldd))
    #identifies final delays which have a single phase delay
    dch_dest2 <- final_delays[which(dd_phases[final_delays] == 1)]
    #creates columns for delayed discharge destinations which have only a single phase
    if (length(dch_dest2) > 0) n <- c(n, paste0(i, "_dd_", dch_dest2))
    #identifies discharge destiantions with multi-phase delays
    dch_dest3 <- final_delays[which(dd_phases[final_delays] > 1)]
    #creates occuancy columns for each delay phase, for discharge destinations with multi-phase delays
    #separate columns for each phase of each unit/service-discharge destination pair
    
    n <- c(n,unlist(sapply(seq(1,length(dch_dest3)),function(y) unlist(sapply(dch_dest3, function(x) paste0(i,"_dd-phase-",y, "_",x))))))
    names(n) <- NULL
    
    #creates columns to count the total number of patients at a unit/service awaiting discharge to a destination with a multi-phase delay
    #(a single column for each unit/service-multi-phase delay destination pair)
    if (length(dch_dest3)>0) n<-c(n,paste0(i,"_dd_", dch_dest3))
    #creates total occupancy columns (active treatment plus delayed transfers/discharges) for each unit/service (one column per unit/service)
    n<-c(n,paste0(i,"_occ"))
    #counts (cumulative) throughput at each unit/service
    n<-c(n,paste0(i,"_tp"))
  }
  
  #records the time (in decimal days) that the system is in the given state (time elapsed until the event of the next row occurs)
  n<- c(n, "time")
  
  ####################### #
  #4.2.2 Inital values ####
  
  #the simulation results matrix "res"
  #starts empty.
  #First row is appended with the initial occupancy, then the matrix builds during the simulation run (one row added per event/iteration of while loop)
  res<-matrix(nrow=0,ncol=length(n))
  colnames(res)<-n
  
  #POPULATE THE FIRST ROW OF RES WITH THE INITIAL OCCUPANCY
  
  #the code section below populates the relevant fields with the initial occupation values
  #it is indexed over the names of res
  #it checks if the name contains "occ" or "act"
  # (if so, this element of res records either a total unit/service occupancy ("occ") or the total of all patients receiving active treatment in the unit ("act"))
  #if it does, return the capacity (bed number) for that unit/service multiplied by the initial occupancy factor (between 0 and 1), rounded up to the nearerst integer
  #then it checks if the name contains the string "_phase"
  #(if so, this element of res records how many patients are currently in a particular phase of active treatment in a unit/service with multi-phase LOS)
  #if it does, check whether the number of phases is a divisor of the capacity for the unit/service multiplied by the initial occupancy factor
  #(by using the modular divisor %% to get the remainder, and checking if the result is zero)
  #if so, then it assigns capacity*initial occupancy/number of phases patients to the phase
  #if not, then
  #if if is the first phase, it assigns the quotient plus the remainder of the capacity*initial occupancy of the unit to this element (%/% + %%)
  #or else if it is a subsequent phase, assign just the quotient of the capaicty*initial occupancy of the unit to this element (just %/%)
  #note that we do not assign any "initial occupancy" of blocked/delayed patients
  #either for transfer between units, or to ultimate discharge destinations (including those with multi-phase delays)
  #it is assumed that at time = 0 the system is operating optimally, with all patients in the system receiving active treamtent in an appropraite unit/service
  #delays and blockages will then accrue from this state during the warm-up period
  #if we do wish to model delays from the outset, we can add extra conditions into the control flow below to split the total occupany between
  #active treatment phases, delayed transfer, and delayed discharge (single-phase or multi-phase)
  #by assigning patients to the columns with "_dd_" and "_dd-phase-" names (and aligning those numbers with the "occ" and "act" counts)
  
  res <- rbind(res, unlist(sapply(1:length(n), function (x){
    if("act" %in% strsplit(n[x],"_")[[1]] || "occ" %in% strsplit(n[x],"_")[[1]]){
      ceiling(beds[[strsplit(n[x],"_")[[1]][1]]]*init_occ)
    } else if (grepl("_phase", n[x]) && ceiling(beds[[strsplit(n[x],"_")[[1]][1]]]*init_occ)%%num_phases[[strsplit(n[x],"_")[[1]][1]]] == 0){
      ceiling(beds[[strsplit(n[x],"_")[[1]][1]]]*init_occ)/num_phases[[strsplit(n[x],"_")[[1]][1]]]
    } else if (grepl("_phase-1", n[x]) && ceiling(beds[[strsplit(n[x],"_")[[1]][1]]]*init_occ)%%num_phases[[strsplit(n[x],"_")[[1]][1]]] != 0){
      ceiling(beds[[strsplit(n[x],"_")[[1]][1]]]*init_occ)%/%num_phases[[strsplit(n[x],"_")[[1]][1]]] +
        ceiling(beds[[strsplit(n[x],"_")[[1]][1]]]*init_occ)%%num_phases[[strsplit(n[x],"_")[[1]][1]]]
    } else if (grepl("_phase", n[x]) && !grepl("_phase-1", n[x])) {ceiling(beds[[strsplit(n[x],"_")[[1]][1]]]*init_occ)%/%num_phases[[strsplit(n[x],"_")[[1]][1]]]
    } else {return(0)}
  }
  )
  )
  )
  
  ######################## #
  #4.3 Event parameters ####
  
  #Named vectors identifying the base parameters (rates) associated with the events which are possible (for the given column names and values in the input CSV)
    #are created here
  #These vectors are then called within the simulation run (while loop), and the parameter values used to calculate event times
    #(which may be dependent on both the base rates and unit occupancy levels at the given iteration of the simulation loop)
  
  #first build vectors of names for arrivals, active service completion and delayed discharge completion
  #and named vectors of their LOS rates, with names corresponding to the names in the event list
  t_ls <- strsplit(n, "_")
  
  #arrivals
  #names for the arrival events
  arrival_names <- paste0("arr_", names(which(!is.na(ar))))
  #take the arrival rates taken directly from the ar vector (defined in the wrapper)
  arrival_rates <- ar[which(!is.na(ar))]
  #assign the arrival event names to the arrival rates vector created above, so it can be referenced in the 
  #event list within the simulation run
  names(arrival_rates) <- arrival_names
  
  #active service completion
  #names for the "complete an active treatment phase" events - multi-phase units
  act_phase_index <- unlist(sapply(1:length(t_ls), function(x) {if(length(t_ls[[x]]) == 2) return(x)}))
  active_phase_names <- unlist(sapply(act_phase_index, function(x) {if (substr(t_ls[[x]][2], 1, 5) == "phase") 
    return(paste0(t_ls[[x]][2], "_", t_ls[[x]][1]))}))
  
  #names for the "complete an active treatment phase" events - single-phase units
  #the test on the number of phases needs to be here - every unit with > 1 phase will  have an "act" column AND "phase-k" columns,
  #but those with only a single phase will ONLY have an "act" column
  one_phase_names <- unlist(sapply(act_phase_index, function(x) {
    if (t_ls[[x]][2] == "act" && num_phases[t_ls[[x]][1]] == 1) 
      return(paste0(t_ls[[x]][2], "_", t_ls[[x]][1]))}))
  
  #create a vector containing the active phase names as they appear in the res matrix
  apn_res_name <- unlist(sapply(act_phase_index, function(x) {if (substr(t_ls[[x]][2], 1, 5) == "phase") 
    return(paste0(t_ls[[x]][1], "_", t_ls[[x]][2]))}))
  #and similarly for the single phase units
  opn_res_name <- unlist(sapply(act_phase_index, function(x) {
    if (t_ls[[x]][2] == "act" && num_phases[t_ls[[x]][1]] == 1) 
      return(paste0(t_ls[[x]][1], "_", t_ls[[x]][2]))}))
  
  #assign the elements of the active_phase_names vector as the names of this res names vector,
  #so that res can be subsetted using active_phase_names to get the current number of patients in that phase, when creating the event list (in the sim loop)
  names(apn_res_name) <- active_phase_names
  names(opn_res_name) <- one_phase_names
  
  #calculate the generic LOS rates for each unit (the same for each phase, since it's Erlang)
  #by dividing the number of phases by the mean LOS (taken from the los vector (defined in the wrapper))
  #Note that within the loop, when the event list is populated, for each unit-phase pair, this figure will be multiplied by the current phase occupancy
  
  active_phase_rates <- unlist(sapply(act_phase_index, function(x) {if (substr(t_ls[[x]][2], 1, 5) == "phase")
    #returns k/los
    return(as.numeric(num_phases[t_ls[[x]][1]])/as.numeric(los[t_ls[[x]][1]]))}))
  
  #assign the active phase completion event names to this vector, so it can be referenced directly in the simulation run
  names(active_phase_rates) <- active_phase_names
  
  #single phase rates
  one_phase_rates <- unlist(sapply(act_phase_index, function(x) {
    if (t_ls[[x]][2] == "act" && num_phases[t_ls[[x]][1]] == 1) 
      return(1/as.numeric(los[t_ls[[x]][1]]))}))
  #assign names to single phase rates vector
  names(one_phase_rates) <- one_phase_names
  
  #delayed discharge completion
  dd_phase_index <- unlist(sapply(1:length(t_ls), function(x) {
    if(length(t_ls[[x]]) == 3 && t_ls[[x]][3] %in% finaldd) return(x)}))
  
  #single phase delayed discharge completion
  #the check that dd_phases for this destination is 1 is necessary because multi-phase destinations also have a "dd" column
  dd_one_phase_names <- unlist(sapply(dd_phase_index, function(x) {
    if (t_ls[[x]][2] == "dd" && dd_phases[[t_ls[[x]][3]]] == 1) return(paste0("ddcomp_", n[[x]] ))}))
  #names as they appear in the res matrix
  ddonepn_res_name <- unlist(sapply(dd_phase_index, function(x) {
    if (t_ls[[x]][2] == "dd" && dd_phases[[t_ls[[x]][3]]] == 1) return(n[[x]] )}))
  #apply names to vector
  names(ddonepn_res_name) <- dd_one_phase_names
  
  #multiple phase delayed discharge
  dd_phase_names <- unlist(sapply(dd_phase_index, function(x) {
    if (substr(t_ls[[x]][2], 1, 8) == "dd-phase") return(paste0("ddcomp_", n[[x]]))}))
  #create a vector containing dd phase names as they appear in the res matrix
  ddpn_res_name <- unlist(sapply(dd_phase_index, function(x) {
    if (substr(t_ls[[x]][2], 1, 8) == "dd-phase") return(n[[x]])}))
  #assign the elements of dd_phase_names as names of this vector so it can be subsetted during event list creation in the sim loop
  names(ddpn_res_name) <- dd_phase_names
  
  #calculate the generic LOS rate for each dd destination (the same for each phase, and regardless of the discharing unit)
  #by dividing the number of delay phases for the destination by the destination's mean delay (taken from the dd vector (defined in the wrapper))
  #Note that within the loop, when the event list is populated, for each unit-phase-destination triple this figure will be multiplied by the 
  #current phase occupancy
  #for single-phase delays
  dd_one_phase_rates <- unlist(sapply(dd_phase_index, function(x) {
    if (t_ls[[x]][2] == "dd" && dd_phases[[t_ls[[x]][3]]] == 1) 
      return(1/as.numeric(dd[t_ls[[x]][3]]))}))
  
  #for multiple-phase delays  
  dd_phase_rates <- unlist(sapply(dd_phase_index, function(x) {if (substr(t_ls[[x]][2], 1, 8) == "dd-phase")
    #returns k/delay
    return(as.numeric(dd_phases[t_ls[[x]][3]])/as.numeric(dd[t_ls[[x]][3]]))}))
  
  #name these vectors so they can be referenced directly in the simulation run
  names(dd_one_phase_rates) <- dd_one_phase_names
  names(dd_phase_rates) <- dd_phase_names
  
  #combine the event names into a single character vector
  event_names <- c(arrival_names, active_phase_names, one_phase_names, dd_one_phase_names, dd_phase_names)
  
  #create an empty event list, the same length as the number of possible (scheduled) events for this simulation scenario
  st <- vector(mode= "list", length = length(event_names))
  #assign the event names to the event list (which can then be referenced in the simulation run loop to pull in the generic rates)
  names(st) <- event_names
  
  
  ###################### #
  #4.4 Simulation run ####
  
  set.seed(SEED)
  #the simulation cycle continues until the (cumulative) time in column "t" exceeds the sum of the input simulation time and warm up time 
  #(in the wrapper)
  while (res[nrow(res),1]<sim_time+wu_time) {
    
    #################### #
    #4.4.1 Event list ####
    
    #4.4.1.1 arrivals ####
    
    #Are we modelling a loss system? If so, we can only have arrivals at units which both take arrivals and are not currently full
    #otherwise we can have an arrival at any unit (but they will be placed in a front door queue when that arrival is processed)
    if(front_door_queue == FALSE) {
      #We are modelling a loss system
      
      #which units can take arrivals (if they're not full)?
      arr_names <- names(ar[which(!is.na(ar))])
      arr_names_occ <- paste0(arr_names,"_occ")
      
      #Which units are not currently full?
      not_full <- arr_names[res[which.max(res[,"t"]), arr_names_occ] < unlist(beds[arr_names])]
      full <- arr_names[res[which.max(res[,"t"]), arr_names_occ] >= unlist(beds[arr_names])]
      #now use this as a conditional statement to decide whether an arrival can happen
      if(length(not_full) > 0){
        #sample the rates for the non-full units
        arrivals <- sapply(arrival_names, function(x) {rexp(1, rate = as.numeric(arrival_rates[[x]]))})
        st[paste0("arr_",not_full)] <- arrivals[paste0("arr_",not_full)]
        st[paste0("arr_",full)] <- NA
      } else {
        st[arrivals] <- NA
      }
    } else {
      #we allow arrivals at any unit which has a defined arrival rate parameter, regardless of its current occupancy
      arrivals <- sapply(arrival_names, function(x) {rexp(1, rate = as.numeric(arrival_rates[[x]]))})
      st[arrival_names] <- arrivals
    }
    
    
    #4.4.1.2 active treatment/service completion ####
    
    #single phases
    onephasecomp <- sapply(one_phase_names, function(x) {rexp(1, rate = res[nrow(res), opn_res_name[x]]*one_phase_rates[x])})
    st[one_phase_names] <- onephasecomp
    
    #multiple phases
    actphasecomp <- sapply(active_phase_names, function(x) {rexp(1, rate = res[nrow(res), apn_res_name[x]]*active_phase_rates[x])})
    st[active_phase_names] <- actphasecomp
    
    #4.4.1.3 delayed discharge completion ####
    
    #single phases
    ddonephasecomp <- sapply(dd_one_phase_names, function(x) {rexp(1, rate = res[nrow(res), ddonepn_res_name[x]]*dd_one_phase_rates[x])})
    st[dd_one_phase_names] <- ddonephasecomp
    
    #multiple phases
    ddphasecomp <- sapply(dd_phase_names, function(x) {rexp(1, rate = res[nrow(res), ddpn_res_name[x]]*dd_phase_rates[x])})
    st[dd_phase_names] <- ddphasecomp
    
    
    ########################### #
    #4.4.2 Action upon event ####
    
    #the first event (smallest number = earliest occurance) is chosen from the list st
    #note that the which.min function returns an index, and it discards NAs and NaNs
    event<-names(which.min(st))
    
    #this adds in the time elapsed until the current event (not the total time, just the time difference), to the "time" column of the previous row of res
    #this is so that total time spent in particular states (e.g. 18 beds occupied at HASU) can be summed in the output measures
    #to get actual time rather than event time measures of occupancy, blockage etc.
    res[nrow(res), "time"] <- min(unlist(st), na.rm = TRUE)
    
    #append the current row (current state) to the end of the matrix res, ready to be modified to refelct the new state following the event
    res <- rbind(res,res[nrow(res),])
    
    #set the event time "t" of the current row of res to the event time of the previous row, plus the time which elapsed between then and the current event
    res[nrow(res),1] <- res[nrow(res)-1,1]+min(unlist(st),na.rm=TRUE)
    
    #
    upstr_tfr<-FALSE
    
    #for the rest of this event (current iteration of the while loop), this value is used to subset for the current row
    current_row <- nrow(res)
    
    #4.4.2.1 arrivals ####
    
    if (substr(event,1,3)=="arr") {
      unit<-unlist(strsplit(event,"_"))
      unit<-unit[length(unit)]
      
      if (res[current_row,which(colnames(res)==paste0(unit,"_occ"))]==beds[[unit]]){
        #unit at capacity - add one to WL
        res[current_row,which(colnames(res)==paste0(unit,"_q"))]<-res[current_row,which(colnames(res)==paste0(unit,"_q"))]+1
      } else if (num_phases[unit] > 1) {
        #multi-phase unit not at capacity - add one to phase 1
        res[current_row,which(colnames(res)==paste0(unit,"_phase-1"))]<-res[current_row,which(colnames(res)==paste0(unit,"_phase-1"))]+1
      } else {
        #single phase unit not at capacity - add one to unit_act
        res[current_row,which(colnames(res)==paste0(unit,"_act"))]<-res[current_row,which(colnames(res)==paste0(unit,"_act"))]+1
      }
    }
    
    #4.4.2.2 active service completion ####

    #IS IT AN ACTIVE TREATMENT PHASE?
    #There are two cases - multi-phase and single phase - and they are named differently in the simulation matrix.
    #Every treatment unit, single or multi-phase, has a field whose name ends with the suffix "act"
    #For single phase units, this corresponds to the number of patients currently in the unit and it is this field which is directly incremented/deincremented
    #For multi-phase units, each treatment phase has a field whose name begins with the prefix "phase" followed by the relevant phase number
    #In this case, it is the fields prefixed "phase" which are directly incremented/deincremented, while the total from all phases (corresponding 
    #to unit active occupancy) is calculated indirectly (along with the total unit occupancy) after all actions resulting from the current event
    #have been processed, and recorded
    
    #4.4.2.2.1 single phase ####
    #(identified by the event name beginning with "act" - n.b. multi-phase active treatment completion event names begin with "phase")
    if (substr(event, 1, 3) == "act"){
      event_substrings <- unlist(strsplit(event,"_"))
      unit <- event_substrings[length(event_substrings)]
      #selects the row from the transition matrix "tr" corresponding to the unit where treatment is being completed
      tmp <- tr[which(rownames(tr)==unit),]#,drop=FALSE]
      tmp <- tmp[, which(!is.na(tmp)==TRUE), drop = FALSE]
      tmp1 <- as.numeric(tmp)
      names(tmp1) <- colnames(tmp)
      #next, partition the unit interval using those probabilities and choose one at random, using a sample from the uniform distribution
      #n.b. the 'fn1' function called here is the one defined at the start of the code, which returns a result as soon as
      #the first input variable (the cumulative sum of trans/dist probs here) is greater than the second input variable (the unif dist sample here)
      dch_dest<-fn1(cumsum(tmp1),runif(1))
      
      #THE LISTS "final" AND "finaldd" are created in the wrapper
      if (dch_dest %in% setdiff(final,finaldd)) {
        #condition true - final discharge no delay: deincrement "act", increment throughput, 
        #and set condition to check for upstream wl/q transfers
        res[current_row,which(colnames(res)==paste0(unit,"_act"))] <-
          res[current_row,which(colnames(res)==paste0(unit,"_act"))]-1
        res[current_row,which(colnames(res)==paste0(unit,"_tp"))] <- 
          res[current_row,which(colnames(res)==paste0(unit,"_tp"))]+1
        upstr_tfr<-TRUE
      } else if (dch_dest %in% finaldd){
        #condition true - delayed discharge, deincrement "act", increment the "waiting at unit for discharge" total column, and
        #for single-phase delay phase discharges, increment the "_dd_" column,
        #for multi-phase delayed discharge destination first delay phase column "_dd-phase-1_"
        res[current_row,which(colnames(res)==paste0(unit,"_act"))]<-
          res[current_row,which(colnames(res)==paste0(unit,"_act"))]-1
        if (dd_phases[[dch_dest]] == 1) {
          res[current_row,which(colnames(res)==paste0(unit,"_dd_",dch_dest))]<-
            res[current_row,which(colnames(res)==paste0(unit,"_dd_",dch_dest))]+1
        } else {
          res[current_row,which(colnames(res)==paste0(unit,"_dd-phase-1_",dch_dest))]<-
            res[current_row,which(colnames(res)==paste0(unit,"_dd-phase-1_",dch_dest))]+1
        }
      } else if (res[current_row,which(colnames(res)==paste0(dch_dest,"_occ"))]==beds[[dch_dest]]){
        #must be a transient state transfer
        #condition true - unit full, deincrement unit "act" and add to relevant delayed transfer ('dc') field
        res[current_row,which(colnames(res)==paste0(unit,"_act"))]<-
          res[current_row,which(colnames(res)==paste0(unit,"_act"))]-1
        res[current_row,which(colnames(res)==paste0(unit,"_dc_",dch_dest))]<-
          res[current_row,which(colnames(res)==paste0(unit,"_dc_",dch_dest))]+1
      } else {
        #condition false - there's space, so so it's a transient state transfer without delay
        #deincrement unit "act", add one to throughput, increment "act" of the new unit (if it's single phase), or "phase-1" (if its multi-phase)
        
        res[current_row,which(colnames(res)==paste0(unit,"_act"))] <-
          res[current_row,which(colnames(res)==paste0(unit,"_act"))]-1
        res[current_row,which(colnames(res)==paste0(unit,"_tp"))] <-
          res[current_row,which(colnames(res)==paste0(unit,"_tp"))]+1
        if (num_phases[[dch_dest]] == 1){
          res[current_row,which(colnames(res)==paste0(dch_dest, "_act"))] <-
            res[current_row,which(colnames(res)==paste0(dch_dest, "_act"))]+1
        } else {
          res[current_row,which(colnames(res)==paste0(dch_dest, "_phase-1"))] <-
            res[current_row,which(colnames(res)==paste0(dch_dest, "_phase-1"))]+1
        }
        #flag that a space has been freed (which queued/delayed patients might be able to move into)
        upstr_tfr<-TRUE
      }
    }
    
    #4.4.2.2.2 multi-phase ####
    if (substr(event, 1, 5) == "phase"){
      #if true, we need to know the name of the unit and the number of the phase so we can carry out appropriate logical tests
      #and indentify which fields of 'res' to increment and deincrement
      event_substrings <- unlist(strsplit(event,"_"))
      unit <- event_substrings[length(event_substrings)]
      phase <- event_substrings[length(event_substrings)-1]
      phase <- unlist(strsplit(phase, "-"))
      phase <- as.numeric(phase[length(phase)])
      
      if (phase < num_phases[[unit]]){
        #condition is true, completion of non-final active treatment phase - shuffle along to the next phase
        res[current_row,which(colnames(res)==paste0(unit,"_phase-", phase))] <- 
          res[current_row,which(colnames(res)==paste0(unit,"_phase-", phase))]-1
        res[current_row,which(colnames(res)==paste0(unit,"_phase-", phase+1))] <- 
          res[current_row,which(colnames(res)==paste0(unit,"_phase-", phase+1))]+1
      } else {
        #condition is false - go to final phase completion
        #if the phase is the final one, we'll need to know where the patient is being transferred/discharged to and store the name as 'dch_dest'
        #first, create a numeric vector of discharge probabilities for the unit
        tmp <- tr[which(rownames(tr)==unit),]
        tmp <- tmp[, which(!is.na(tmp)==TRUE), drop = FALSE]
        tmp1 <- as.numeric(tmp)
        names(tmp1) <- colnames(tmp)
        #next, partition the unit interval using those probabilities and choose one at random, using a sample from the uniform distribution
        #n.b. the 'fn1' function called here is the one defined at the start of the code, which returns a result as
        #the first input variable (the cumulative sum of trans/dist probs here) is greater than the second input variable (the uniform distribution sample here)
        dch_dest<-fn1(cumsum(tmp1),runif(1))
        
        #the lists "final" and "finaldd" are defined in the wrapper
        if (dch_dest %in% setdiff(final,finaldd)) {
          #condition true - final discharge no delay: deincrement final phase, increment throughput, 
          #and set condition to check for upstream wl/q transfers
          
          res[current_row,which(colnames(res)==paste0(unit,"_phase-", phase))] <-
            res[current_row,which(colnames(res)==paste0(unit,"_phase-", phase))]-1
          res[current_row,which(colnames(res)==paste0(unit,"_tp"))] <- 
            res[current_row,which(colnames(res)==paste0(unit,"_tp"))]+1
          #flag that a space has been freed which a patient waiting upstream might move into
          upstr_tfr<-TRUE
        } else if (dch_dest %in% finaldd){
          #condition true - delayed discharge, deincrement final phase, increment either "_dd_" or "_dd-phase-1_" of the discharge destination
          
          res[current_row,which(colnames(res)==paste0(unit,"_phase-", phase))]<-
            res[current_row,which(colnames(res)==paste0(unit,"_phase-", phase))]-1
          if (dd_phases[[dch_dest]] == 1) {
            #single phase delay for discharge destination
            #just increment the _dd_ column
            res[current_row,which(colnames(res)==paste0(unit,"_dd_",dch_dest))]<-
              res[current_row,which(colnames(res)==paste0(unit,"_dd_",dch_dest))]+1
          } else {
            #multi-phase delay for discharge destination
            #increment the first delay phase column
            
            res[current_row,which(colnames(res)==paste0(unit,"_dd-phase-1_",dch_dest))]<-
              res[current_row,which(colnames(res)==paste0(unit,"_dd-phase-1_",dch_dest))]+1
          }
        } else if (res[current_row,which(colnames(res)==paste0(dch_dest,"_occ"))]==beds[[dch_dest]]){
          #must be a transient state transfer
          #condition true - unit full, deincrement final phase and add to relevant delayed transfer ('dc') field
          
          res[current_row,which(colnames(res)==paste0(unit,"_phase-", phase))]<-
            res[current_row,which(colnames(res)==paste0(unit,"_phase-", phase))]-1
          res[current_row,which(colnames(res)==paste0(unit,"_dc_",dch_dest))]<-
            res[current_row,which(colnames(res)==paste0(unit,"_dc_",dch_dest))]+1
        } else {
          #condition false - there's space, so so it's a transient state transfer without delay
          #deincrement final phase, add one to throughput, increment phase 1 of the new unit,
          #and set condition to check for upstream wl/q transfers
          
          res[current_row,which(colnames(res)==paste0(unit,"_phase-", phase))] <-
            res[current_row,which(colnames(res)==paste0(unit,"_phase-", phase))]-1
          res[current_row,which(colnames(res)==paste0(unit,"_tp"))] <-
            res[current_row,which(colnames(res)==paste0(unit,"_tp"))]+1
          if (num_phases[[dch_dest]] == 1){
            res[current_row,which(colnames(res)==paste0(dch_dest, "_act"))] <-
              res[current_row,which(colnames(res)==paste0(dch_dest, "_act"))]+1
          } else {
            
            res[current_row,which(colnames(res)==paste0(dch_dest, "_phase-1"))] <-
              res[current_row,which(colnames(res)==paste0(dch_dest, "_phase-1"))]+1
          }
          #flag that a space has been freed, which a patient waiting upstream might move into
          upstr_tfr<-TRUE
        }
      }
    }
    
    ######################################## #
    #4.4.2.3 Delayed discharge completion ####
    
    if (substr(event,1,6)=="ddcomp") {
      #if true, then we need to know the name of the discharging unit, the name of the discharge destination, and the number of delay phases
      #we should then proceed with a simplified version of the active treatment phase completion code
      #(with just two options - shuffle along (deincrement current phase, increment next), 
      #or actually leave(deincrement current phase, increment unit throughput, and set upstream transfer considerations to true))
      event_components <- unlist(strsplit(event,"_"))
      unit <- event_components[length(event_components)-2]
      dch_dest <- event_components[length(event_components)]
      
      #the code above has identified that the event is a delayed discharge phase completion and
      #saved the name of the discharging unit as the character 'unit', and
      #saved the name of the discharge destination as the character 'dch_dest'
      
      #single phase
      if (dd_phases[[dch_dest]] == 1){
        #it's a one-phase discharge, so we just deincrement the waiting list and increment the throughput at the discharing unit, then check for upstream transfers
        res[current_row,which(colnames(res)==paste0(unit,"_dd_", dch_dest))] <-
          res[current_row,which(colnames(res)==paste0(unit,"_dd_", dch_dest))] -1
        res[current_row,which(colnames(res)==paste0(unit,"_tp"))] <- 
          res[current_row,which(colnames(res)==paste0(unit,"_tp"))] +1
        #flag that a space has been freed, which a patient waiting upstream might move into
        upstr_tfr<-TRUE
      } else {
        #multi phase
        phase_name <- event_components[length(event_components)-1]
        phase_components <- unlist(strsplit(phase_name, "-"))
        phase <- as.numeric(phase_components[length(phase_components)])
        #the above code has saved the number of the current delay phase as the number 'phase'
        #now we need to check whether it's the final phase or not, and either shuffle along or complete the discharge
        
        if (phase < dd_phases[[dch_dest]]){
          #condition true - it's an intermediary phase, so shuffle along
          res[current_row, which(colnames(res) == paste0(unit, "_dd-phase-", phase, "_", dch_dest))] <- 
            res[current_row, which(colnames(res) == paste0(unit, "_dd-phase-", phase, "_", dch_dest))] - 1
          res[current_row, which(colnames(res) == paste0(unit, "_dd-phase-", phase+1, "_", dch_dest))] <- 
            res[current_row, which(colnames(res) == paste0(unit, "_dd-phase-", phase+1, "_", dch_dest))] +1
        } else {
          #condition false - final phase, so complete discharge
          
          res[current_row,which(colnames(res)==paste0(unit,"_dd-phase-", phase, "_", dch_dest))] <-
            res[current_row,which(colnames(res)==paste0(unit,"_dd-phase-", phase, "_", dch_dest))] -1
          res[current_row,which(colnames(res)==paste0(unit,"_tp"))] <- 
            res[current_row,which(colnames(res)==paste0(unit,"_tp"))] +1
          #flag that a space has been freed which a patient waiting upstream might move into
          upstr_tfr<-TRUE
        }
      }
    }
    
    ########################### #
    #4.4.3 Upstream transers ####
    
    #Has the event just completed freed a space which a patient who is in a queue, or is delayed awaiting transfer, can move into?
    #e.g. if we have just had an active service completion or a delayed discharge completion, has this freed a space in a
    #unit which was full? If so, are there any patients waiting to be admitted/transferred to that unit? If there are, shuffle them along.
    
    if (upstr_tfr==TRUE) {
      #condition true - the most recent event has resulted in a space at that unit being freed
      #now we need to check whether there are any patients waiting to occupy it, and if there are shuffle them down into it
      t_ls<-strsplit(colnames(res),"_")
      loop<-TRUE
      
      while(loop==TRUE) {
        #first time round, this loop looks for patients waiting for the place freed up by the scheduled event and if there are any, shuffles them down
        #if there's no-one waiting, no shuffling happens
        #if there is, then a place has been freed at the unit the patient was shuffled from
        #so the whole thing happens again, shifted upstream by one unit
        #as soon as the loop fails to find anyone waiting at the next step upstream from the most recently freed space, it sets "loop==FALSE" and stops
        
        #find all referring units OR WL
        ref_ind <- unlist(sapply(1:length(t_ls),function(x) if(t_ls[[x]][length(t_ls[[x]])]==unit & length(intersect(t_ls[[x]],"dc"))>0) return(x)))
        ref_unit_ls <- res[current_row,ref_ind]
        if (paste0(unit,"_q") %in% colnames(res)) ref_unit_ls<-c(ref_unit_ls,res[current_row,which(colnames(res)==paste0(unit,"_q"))])
        #choose the name of the unit with the largest waiting list (if there's a tie, the first unit listed in the vector will be chosen)
        #return the name of the waiting list, e.g. "asu_dc_rhb" or "hasu_q"
        ref_unit <- names(which(ref_unit_ls==max(ref_unit_ls))[1])
        
        
        if (sum(ref_unit_ls) == 0){
          #condition true - there are no patients waiting upstream for the space just freed, end the loop
          loop <- FALSE
        } else if (ref_unit == paste0(unit, "_q")) {
          #patients are waiting for the newly freed space
          #condition true - admission from a front door queue (no upstream space is freed here so the loop halts after admission)
          if (num_phases[[unit]] == 1){
            #admission to a single phase unit
            res[current_row,which(colnames(res)==paste0(unit,"_act"))] <- 
              res[current_row,which(colnames(res)==paste0(unit,"_act"))]+1
            res[current_row,which(colnames(res)==paste0(unit,"_q"))] <- 
              res[current_row,which(colnames(res)==paste0(unit,"_q"))]-1
          } else {
            #admission to first phase of a multi phase unit
            
            res[current_row,which(colnames(res)==paste0(unit,"_phase-1"))] <- 
              res[current_row,which(colnames(res)==paste0(unit,"_phase-1"))]+1
            res[current_row,which(colnames(res)==paste0(unit,"_q"))] <- 
              res[current_row,which(colnames(res)==paste0(unit,"_q"))]-1
          }
          #either way, exit the loop
          loop<-FALSE
        } else {
          #patients are waiting for the new freed space, and since the previous condition was false, they must be blocked at a transient unit 
          #deincrement referring unit dc list, increment referring unit tp
          ref_unit<-unlist(strsplit(ref_unit,"_"))[1]
          if (num_phases[[unit]] == 1){
            #admission to a single phase unit
            res[current_row,which(colnames(res)==paste0(unit,"_act"))] <- 
              res[current_row,which(colnames(res)==paste0(unit,"_act"))]+1
            res[current_row,which(colnames(res)==paste0(ref_unit,"_dc_",unit))] <- 
              res[current_row,which(colnames(res)==paste0(ref_unit,"_dc_",unit))]-1
            res[current_row,which(colnames(res)==paste0(ref_unit,"_tp"))] <- 
              res[current_row,which(colnames(res)==paste0(ref_unit,"_tp"))]+1
            
          } else {
            #admission to a multi phase unit
            res[current_row,which(colnames(res)==paste0(unit,"_phase-1"))] <- 
              res[current_row,which(colnames(res)==paste0(unit,"_phase-1"))]+1
            #CHANGE 20171205 - MOVED THE "_dc_" DEINCREMENT AND THE "tp" INCREMENTS HERE FROM BELOW
            res[current_row,which(colnames(res)==paste0(ref_unit,"_dc_",unit))] <- 
              res[current_row,which(colnames(res)==paste0(ref_unit,"_dc_",unit))]-1
            res[current_row,which(colnames(res)==paste0(ref_unit,"_tp"))] <- 
              res[current_row,which(colnames(res)==paste0(ref_unit,"_tp"))]+1
          }
          
          #this has freed a bed at the referring unit
          #so we move upstream by setting the unit which has just shuffled down a patient (the ref_unit) to now be the current unit (unit)
          #and go through the loop again, to check if any patients are waiting for that space
          #and continue to move upstream until the chain is completed, and there are no more patients waiting to shuffle
          #into the spaces freed by the current event list event
          unit<-ref_unit
        }
      }
    }
    
    ########################## #
    #4.4.4 Occupancy counts ####
    
    #All events in this iteration of the simulation run (triggered by the most recent event drawn from the event list) have been completed
    #counts of total occupancy ("occ"), total number of patients receiving active treatment ("act"), and total numbers of patients delayed awaiting discharge ("dd")
      #at each unit are now calculated before moving on to the next event
    
    t_ls <- strsplit(colnames(res), "_")
    
    #generates totals for occupancy, active treatment and delayed discharge (by destination) for each unit at the current state
    for (i in trans){
      #create an index of active treatment phase fields for unit i
      tmp_act_phases <- unlist(sapply(1:length(t_ls), 
                                      function(x) {
                                        if(t_ls[[x]][1] == i && length(t_ls[[x]]) > 1 && substring(t_ls[[x]][2],1,5)=="phase") return(x)}))
      
      #create an index of all delayed transfer to transient state fields from unit i (at most one for each unit pair - no phases here)
      tmp_dc <- unlist(sapply(1:length(t_ls), 
                              function(x) {
                                if(t_ls[[x]][1] == i && length(t_ls[[x]]) > 1 && t_ls[[x]][2] == "dc") return(x)}))
      
      #create an index of all delayed discharge fields from unit i
      #(this index does not discriminate between discharge destinations - it will be subsetted in the nested loop to split by destination)
      #for multi phase delayed discharges
      tmp_dd_phases <- unlist(sapply(1:length(t_ls), 
                                     function(x) {
                                       if(t_ls[[x]][1] == i && length(t_ls[[x]]) > 1 && substring(t_ls[[x]][2],1,8)=="dd-phase") return(x)}))
      #for single phase delayed discharges
      
      tmp_dd_one_phase_index <- unlist(sapply(1:length(t_ls), function(x) {if(length(t_ls[[x]]) == 3) return(x)}))
      tmp_dd_one_phase_index <- unlist(sapply(tmp_dd_one_phase_index, function(x) {if(t_ls[[x]][2] == "dd") return(x)}))
      
      tmp_dd_one_phase <- unlist(sapply(tmp_dd_one_phase_index, 
                                        function(x) {
                                          if(t_ls[[x]][1] == i && dd_phases[[t_ls[[x]][3]]] == 1) return(x)}))
      
      #for single phase delayed discharges, the count by destination has already been generated above
      #for multi phase discharges, sum the phases to get totals of delayed discharges from i split by destination
      for (j in finaldd){
        if (dd_phases[[j]] > 1){
          #condition false - it's a single phase discharge to leave it alone
          #condition true - create an index of delayed discharge destination fields from i, 
          #by subsetting the event list using the index of all delayed discharges from unit i
          tmp_dd_dest_phases <- unlist(sapply(tmp_dd_phases, 
                                              function(x) {
                                                if(length(t_ls[[x]]) > 1 && t_ls[[x]][length(t_ls[[x]])] == j) return(x)}))
          #sum delayed discharge phases from unit i to destination j to get the delay from i to j total (at this state)
          if (length(tmp_dd_dest_phases)>0){
            res[current_row, paste0(i, "_dd_", j)] <- sum(res[current_row, tmp_dd_dest_phases])
          }
        }
      }
      
      #sum active treatment phases for i to get the active treamtent total
      #if there's only one phase, it has already been populated in the run above, so leave it alone
      if (num_phases[[i]] > 1){
        res[current_row, paste0(i, "_act")] <- sum(res[current_row, tmp_act_phases])
      }
      
      #sum active the treamtent total, delayed transfers, and delayed discharges to all destinations from i, to get the current total occupancy at i 
      res[current_row, paste0(i, "_occ")] <- res[current_row, paste0(i, "_act")] + sum(res[current_row, tmp_dc]) + 
        sum(res[current_row, tmp_dd_phases]) + sum(res[current_row, tmp_dd_one_phase])
      
      #generates totals of delayed discharges from i split by destination
      for (j in finaldd){
        #create an index of delayed discharge destination fields from i by subsetting the event list using the index of 
        #all delayed discharges from unit i
        tmp_dd_dest_phases <- unlist(sapply(tmp_dd_phases, 
                                            function(x) {if(length(t_ls[[x]]) > 1 & t_ls[[x]][length(t_ls[[x]])] == j) return(x)}))
        
        #sum delayed discharge phases from unit i to destination j to get the delay from i to j total (at this state)
        if (length(tmp_dd_dest_phases)>0){
          res[current_row, paste0(i, "_dd_", j)] <- sum(res[current_row, tmp_dd_dest_phases])
        }
      }
    }
    
    ################ #
    #4.4.5 Checks ####
    
    #check occupancy doesn't exceed capacity
    for (i in trans) {
      if (res[current_row,which(colnames(res)==paste0(i,"_occ"))]>beds[[i]]) stop("Occupancy exceeds capacity")
    }
    #check if occupancy is non integer
    for (i in trans) {
      if (res[current_row,which(colnames(res)==paste0(i,"_occ"))]!=round(res[current_row,which(colnames(res)==paste0(i,"_occ"))])) stop("Occupancy is non-integer")
    }
    #check if occupancy is negative
    for (i in trans){
      if (res[current_row, which(colnames(res) == paste0(i, "_occ"))] < 0) stop("Occupancy is negative")
    }
    #check if active treatment count is negative
    for (i in trans){
      if (res[current_row, which(colnames(res) == paste0(i, "_act"))] < 0) stop("Number of patients in active treatment is negative")
    }
    
    #check there is not spare capacity and a WL for units at same time
    for (i in trans) {
      bedsfree<-beds[[i]]-res[current_row,which(colnames(res)==paste0(i,"_occ"))]
      wls_ind<-unlist(sapply(1:length(t_ls),function(x) {
        if(t_ls[[x]][length(t_ls[[x]])]==i && length(intersect(t_ls[[x]],"dc"))>0) return(x)}))
      wls<-sum(res[current_row,wls_ind])
      if (bedsfree>0 & wls>0) stop("spare capacity and WL")
    }
    
  }
  
  ####################################### #
  #4.5 Outputs (single simulation run) ####
  
  #4.5.1 Simulation matrix ####
  res <- res[which(!is.na(res[,"t"])),]
  write.csv(res,file=paste0(getwd(),"/",folder_name,"/res_",runs,".csv"))
  
  #4.5.2 Plots ####
  
  #4.5.2.1 Plot 1 (occupancy and blocked) ####
  pdf(paste0(getwd(),"/",folder_name,"/p1_",runs,".pdf"),width=8,height=8)
  states<-1+1:length(c(trans,finaldd))
  names(states)<-c(trans,finaldd)
  tmp1<-res[,1]-wu_time
  for (i in trans) {
    tmp2<-res[,which(colnames(res)==paste0(i,"_occ"))]
    #Q for ....
    tmp5<-unlist(sapply(1:length(t_ls),function(x) if(t_ls[[x]][1]==i & length(intersect(t_ls[[x]],"dd"))>0) return(x)))
    tmp5<-res[,tmp5,drop=FALSE]
    tmp6<-unlist(sapply(1:length(t_ls),function(x) if(t_ls[[x]][1]==i & length(intersect(t_ls[[x]],"dc"))>0) return(x)))
    tmp6<-res[,tmp6,drop=FALSE]
    tmp6a<-rowSums(cbind(tmp5,tmp6))
    #plot
    plot(tmp1,tmp2,type="l",col="black",ylim=c(0,1.2*max(tmp2,beds[[i]])),xlab="Time (days)"
         ,ylab="Number of patients / referrals",lwd=2,main=i)
    abline(h=beds[[i]],col="darkgrey",lty="dashed")
    abline(v=0,col="darkgrey")
    lines(tmp1,tmp6a,col="red")
    legend("topleft",c("Occupied beds","Bed capacity","Beds blocked"),col=c("black","darkgrey","red"),lwd=c(2,rep(1,2))
           ,lty=c("solid","dashed","solid"),bg="white")
  }
  dev.off()
  
  #4.5.2.2 Plot 2 (occupancy and queues) ####
  pdf(paste0(getwd(),"/",folder_name,"/p2_",runs,".pdf"),width=8,height=8)
  states<-1+1:length(c(trans,finaldd))
  names(states)<-c(trans,finaldd)
  tmp1<-res[,1]-wu_time
  for (i in trans) {
    tmp2<-res[,which(colnames(res)==paste0(i,"_occ"))]
    #Q for ....
    tmp5<-unlist(sapply(1:length(t_ls),function(x) if(t_ls[[x]][1]==i & length(intersect(t_ls[[x]],"dd"))>0) return(x)))
    tmp5<-res[,tmp5,drop=FALSE]
    tmp6<-unlist(sapply(1:length(t_ls),function(x) if(t_ls[[x]][1]==i & length(intersect(t_ls[[x]],"dc"))>0) return(x)))
    tmp6<-res[,tmp6,drop=FALSE]
    #Q at ....
    tmp7<-unlist(sapply(1:length(t_ls),function(x) if(t_ls[[x]][length(t_ls[[x]])]==i & length(intersect(t_ls[[x]],"dc"))>0) return(x)))
    tmp7<-res[,tmp7,drop=FALSE]
    #plot
    plot(tmp1,tmp2,type="l",col="black",ylim=c(0,1.2*max(tmp2,beds[[i]],tmp5,tmp6,tmp7)),xlab="Time (days)"
         ,ylab="Number of patients / referrals",lwd=2,main=i,col.main=states[which(names(states)==i)])
    abline(h=beds[[i]],col="darkgrey",lty="dashed")
    abline(v=0,col="darkgrey")
    tmp5a<-NULL
    tmp6a<-NULL
    tmp7a<-NULL
    if (ncol(tmp5)>0) {
      tmp5a<-sapply(1:ncol(tmp5),function(x) states[which(names(states)==strsplit(colnames(tmp5)[x],"_")[[1]][length(strsplit(colnames(tmp5)[x],"_")[[1]])])])
      for (j in 1:ncol(tmp5)) lines(tmp1,tmp5[,j],col=tmp5a[j],lty="dotted")
    }
    if (ncol(tmp6)>0) {
      tmp6a<-sapply(1:ncol(tmp6),function(x) states[which(names(states)==strsplit(colnames(tmp6)[x],"_")[[1]][length(strsplit(colnames(tmp6)[x],"_")[[1]])])])
      for (j in 1:ncol(tmp6)) lines(tmp1,tmp6[,j],col=tmp6a[j])
    }
    if (ncol(tmp7)>0) {
      tmp7a<-sapply(1:ncol(tmp7),function(x) states[which(names(states)==strsplit(colnames(tmp7)[x],"_")[[1]][1])])
      for (j in 1:ncol(tmp7)) lines(tmp1,tmp7[,j],col=tmp7a[j])
    }
    leg_Qto_dd<-if (ncol(tmp5)>0) paste("Q for",sapply(strsplit(colnames(tmp5),"_"),function(x) x[length(x)])) else NULL
    leg_Qto_dc<-if (ncol(tmp6)>0) paste("Q for",sapply(strsplit(colnames(tmp6),"_"),function(x) x[length(x)])) else NULL
    leg_Qat_dc<-if (ncol(tmp7)>0) paste("Q at",sapply(strsplit(colnames(tmp7),"_"),function(x) x[1])) else NULL
    leg<-c(leg_Qto_dd,leg_Qto_dc,leg_Qat_dc)
    tcol<-as.numeric(c(tmp5a,tmp6a,tmp7a))
    legend("topleft",c("Occupied beds","Bed capacity",leg),col=c("black","darkgrey",tcol),lwd=c(2,rep(1,1+length(leg)))
           ,lty=c("solid","dashed",rep("dotted",length(tmp5a)),rep("solid",length(leg)-length(tmp5a))),bg="white")
  }
  dev.off()
  
  
  #4.5.3 Performance results (summary statistics) ####
  
  tres<-res[which(res[,1]>=wu_time),]
  tres[,1]<-tres[,1]-wu_time
  
  RES<-matrix(nrow=0,ncol=15)
  for (i in trans) {
    tRES<-rep(NA,13)
    
    #4.5.3.1 Time at occupancy/blockage level calculations ####
    
    #create a temporary dataframe ("tmpocctimes") to store the total amount of time spent at each occupancy level,
    #and the cumulative amount of time ("cu_time") spent at or below that occupancy level
    #for each unit which has an occupancy
    #the entries of this dataframe are used in the occupancy measures below
    
    tempocc1 <- tres[, c(paste0(i, "_occ"), "time")]
    tempocc2 <- tapply(tempocc1[, "time"], tempocc1[, paste0(i, "_occ")], sum)
    tmpocctimes <- data.frame(occ = numeric(length(tempocc2)), time = numeric(length(tempocc2)), cu_time = numeric(length(tempocc2)))
    tmpocctimes[, "occ"] <- as.numeric(names(tempocc2))
    tmpocctimes[, "time"] <- tempocc2
    tmpocctimes$cu_time <- cumsum(tmpocctimes$time)
    
    #create a temporary dataframe ("tempblocktimes") to store the total amount of time that the given number of beds/service channels are blocked
    #and the cumulative amount of time ("cu_time") spent at or below that blockage level
    #for each unit which has an occupancy
    #the entries of this dataframe are used in the "beds blocked" measures below
    
    tempblock1 <- tres[, c(paste0(i, "_occ"), paste0(i, "_act"), "time")]
    blocked <- tempblock1[, paste0(i, "_occ")] - tempblock1[, paste0(i, "_act")]
    tempblock1 <- cbind(tempblock1, blocked)
    tempblock2 <- tapply(tempblock1[, "time"], tempblock1[, "blocked"], sum)
    tempblocktimes <- data.frame(block = numeric(length(tempblock2)), time = numeric(length(tempblock2)), cu_time = numeric(length(tempblock2)))
    tempblocktimes[, "block"] <- as.numeric(names(tempblock2))
    tempblocktimes[, "time"] <- tempblock2
    tempblocktimes$cu_time <- cumsum(tempblocktimes$time)
    
    #create a temporary dataframe ("tmpqtimes") to store the total amount of time the waiting list is at each size
    #and the cumulative amount of time spent with a waiting list that size or smaller
    #for each unit which can have a waiting list (either at arrival, or at another unit)
    #the entries in this dataframe are used in the waiting list measures below
    
    #creates an index of the columns which count delayed transfers into this unit from upstream units
    tmpdelayqs <- unlist(sapply(1:length(t_ls),function(x) if(t_ls[[x]][length(t_ls[[x]])]==i & length(intersect(t_ls[[x]],"dc"))>0) return(x)))
    #if the unit can have both a front door queue and queues at other units, sum all of them to get the total queue 
    if (
      length(intersect(paste0(i, "_q"), colnames(tres))) > 0 && length(tmpdelayqs) > 0) {
      tempq1 <- tres[, c(grep(paste0(i, "_q"), colnames(tres)), tmpdelayqs, grep("time", colnames(tres)))]
      total_q <- tempq1[,paste0(i, "_q")] + tempq1[, n[tmpdelayqs]]
      tempq2 <- tapply(tempq1[, "time"], total_q, sum)
      #if the unit can have a front door queue, but not queues at other units, return just the front door queue as the total queue
    } else if (length(intersect(paste0(i, "_q"), colnames(tres))) > 0) {
      tempq1 <- tres[, c(paste0(i, "_q"), "time")]
      total_q <- tempq1[,paste0(i, "_q")]
      tempq2 <- tapply(tempq1[, "time"], total_q, sum)
    } else {
      #if the unit can have queues at other units, but not a front door queue, sum just the queues at other units as the total queue
      tempq1 <- tres[, c(tmpdelayqs, grep("time", colnames(tres)))]
      total_q <- tempq1[, 1]
      tempq2 <- tapply(tempq1[, "time"], total_q, sum)
    }
    
    tempqtimes <- data.frame(queue = numeric(length(tempq2)), time = numeric(length(tempq2)), 
                             cu_time = numeric(length(tempq2)))
    tempqtimes[, "queue"] <- as.numeric(names(tempq2))
    tempqtimes[, "time"] <- tempq2
    tempqtimes$cu_time <- cumsum(tempqtimes$time)
    
    
    #4.5.3.2 Summary statistics matrix ####
    

    # i. mean occupancy # beds
    #proportion of time at each occupancy assigned to vector temp4
    temp4 <- tmpocctimes[, "time"]/tres[nrow(tres), 1]
    #calculate mean occ by summing product of occupancy and temp4
    tRES[1]<- sum(tmpocctimes[, "occ"]*temp4)
    
    # ii. mean occupancy % beds of capacity
    #calculate mean occ as % of capacity by summing product of proportion of capacity for each occupancy with proportion fo time spent at occupancy
    tRES[2]<-sum((tmpocctimes[, "occ"]/beds[[i]])*temp4)
    
    # iii. 10 quantile occupancy
    tRES[3] <- min(tmpocctimes[["occ"]][which(tmpocctimes[["cu_time"]] >= 0.10*tres[nrow(tres), 1])])
    
    # iv. 50 quantile occupancy 
    tmp <- tres[, which(colnames(tres) == paste0(i, "_occ"))]
    tRES[4] <- min(tmpocctimes[["occ"]][which(tmpocctimes[["cu_time"]] >= 0.50*tres[nrow(tres), 1])])
    
    # v. 85 quantile occupancy
    tmp <- tres[, which(colnames(tres) == paste0(i, "_occ"))]
    tRES[5] <- min(tmpocctimes[["occ"]][which(tmpocctimes[["cu_time"]] >= 0.85*tres[nrow(tres), 1])])
    
    # vi. 90 quantile occupancy
    tmp <- tres[, which(colnames(tres) == paste0(i, "_occ"))]
    tRES[6] <- min(tmpocctimes[["occ"]][which(tmpocctimes[["cu_time"]] >= 0.90*tres[nrow(tres), 1])])
    
    # vii. 95 quantile occupancy  
    tmp <- tres[, which(colnames(tres) == paste0(i, "_occ"))]
    tRES[7] <- min(tmpocctimes[["occ"]][which(tmpocctimes[["cu_time"]] >= 0.95*tres[nrow(tres), 1])])
    
    # viii. 99 quantile occupancy
    tmp <- tres[, which(colnames(tres) == paste0(i, "_occ"))]
    tRES[8] <- min(tmpocctimes[["occ"]][which(tmpocctimes[["cu_time"]] >= 0.99*tres[nrow(tres), 1])])
    
    # ix. % time at full occupancy
    #calculate % time at full occ directly as amount of time at full occ divided by total time
    if (length(tmpocctimes[which(tmpocctimes[, "occ"] == beds[[i]]), "occ"] > 0 )) {
      tRES[9] <- tmpocctimes[which(tmpocctimes[, "occ"]/beds[[i]] ==1), "time"]/tres[nrow(tres), 1]
    } else { 
      tRES[9] <- 0}
    
    # x. mean beds blocked
    tempblkprop <- tempblocktimes[, "time"]/tres[nrow(tres), 1]
    tRES[10]<- sum(tempblocktimes[, "block"]*tempblkprop)
    
    # xi. 95 quant beds blocked
    tRES[11]<- min(tempblocktimes[["block"]][which(tempblocktimes[["cu_time"]] >= 0.95*tres[nrow(tres), 1])])
    
    # xii. mean % beds blocked
    tRES[12] <- sum((tempblocktimes[, "block"]/beds[[i]])*tempblkprop)
    
    # xiii. 95 quant % beds blocked
    tRES[13] <- min(tempblocktimes[["block"]][which(tempblocktimes[["cu_time"]] >= 0.95*tres[nrow(tres), 1])])/beds[[i]]
    
    # save to file
    RES<-rbind(RES,c(runs,i,tRES))
    
  }
  
  return(RES)
}

#################################### #
#5. PARALLELISATION AND ITERATION ####
#################################### #

#In this code section, the simulation function defined above is called and run the number of times specified in the "nruns" input
#Individual runs are carried out in parallel across different CPU cores (on a single Windows computer) then recombined before summary measures are calculated

start.time<-Sys.time()

cl<-makeCluster(7)
registerDoSNOW(cl)
RESULTS<-foreach(runs=1:nruns,.combine="rbind") %dopar% {
  SEED<-seedstart+runs
  RES<-sim_run()
  print(paste0("Results for run ",runs))
  return(RES)
}
stopCluster(cl)

end.time<-Sys.time()
time.taken<-end.time-start.time
print(time.taken)

#################################################### #
#6. AGGREGATION OF RESULTS ACROSS SIMULATION RUNS ####
#################################################### #

#This section calculates summary statistics (occupancy quantiles, percentage of time at full capacity, blockages, means) aggregated across all of the runs
  #in the simulation, and writes them to CSV
#note that results for each indivudal run are written to CSV within the function above

n<-c("run","unit","occ_mean","occ_mean_%_cap", "occ_10", "occ_50", "occ_85", "occ_90", "occ_95",
     "occ_99", "full_cap_%", "blocked_mean","blocked_95","blocked_mean_%_occ","blocked_95_%_occ")
colnames(RESULTS)<-n
write.csv(RESULTS,file=paste0(getwd(),"/",folder_name,"/ALL_RES.csv"))

sumRESULTS<-sapply(1:length(trans),function(x) {
  tmp<-RESULTS[which(RESULTS[,2]==trans[x]),3:15]
  class(tmp)<-"numeric"
  colMeans(tmp)
})
sumRESULTS<-t(sumRESULTS)
sumRESULTS<-cbind(trans,sumRESULTS)
colnames(sumRESULTS)<-n[-1]
#write.csv(sumRESULTS,file=paste0(getwd(),"/",folder_name,"/sum_ALL_RES.csv"))

TsumRESULTS<-as.data.frame(sumRESULTS[,c(1,10,2,6,8,9,13,14)],stringsAsFactors=FALSE)
TsumRESULTS$occ_mean<-round(as.numeric(TsumRESULTS$occ_mean),1)
TsumRESULTS$occ_85<-round(as.numeric(TsumRESULTS$occ_85),1)
TsumRESULTS$occ_95<-round(as.numeric(TsumRESULTS$occ_95),1)
TsumRESULTS$occ_99<-round(as.numeric(TsumRESULTS$occ_99),1)
TsumRESULTS$"full_cap_%"<-round(as.numeric(TsumRESULTS$"full_cap_%"),2)
TsumRESULTS$"blocked_mean_%_occ"<-round(as.numeric(TsumRESULTS$"blocked_mean_%_occ"),2)
TsumRESULTS$"blocked_95_%_occ"<-round(as.numeric(TsumRESULTS$"blocked_95_%_occ"),2)
write.csv(TsumRESULTS,file=paste0(getwd(),"/",folder_name,"/sum_ALL_RES.csv"))

