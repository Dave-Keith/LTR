# This function is the workhorse to make the model 'fit' to the observed time series, to do so we adjust the fecundity or survivorship variables iteratively until
# the observed population abundance is within X% of the VPA model.

# All the data you need to make it dance...
# years:          The years you are running the backwards calculation for
# step.size:      How big is the difference (proportional) in the simulation step.  Default = 0.05 which is about a 5% change for each step
# tuner:          What are we tuning, options are 'fecund' for fecundity, 'fec_nm' for fecundity and natural mortality, 'z' for natural mortality and 
#                 fishing mortality, 'nm' for natural mortality, 
# direction       Allow for forwards and backwards simulations
# ages:           The age classes of your population
# nm:             Natural mortality, this wants the instantaneous natural mortality, not the proportional. Should be a matrix the rows are different years
#                 and the columns are the different ages.
# fecund:         Fecundity is the number of recruits produced by the average individual in the age classes.  Should be a matrix the rows are different years
#                 and the columns are the different ages.
# fm:             Fishing mortality, should be a matrix with the columns being the age classes and rows being years. Should be instantaneous fishing mortality
# N.init          The abundance in year 1
# abund.ts:       The abundance time series from the assessment
# abund.age       The abundance by age, should be a matrix with rows being years and columns being age classes
# catch.age       The catch by age, should be a matrix with rows being years and columns being age classes. This is in numbers not biomass
# preload:        if you have already loaded the functions you need...

fast.tunes<-function(years,step.size=0.05,tuner="m",  direction = 'backwards',ages = NULL,
                   nm=NULL,fecund = NULL,fm = NULL,N.init = NULL,abund.ts = NULL,abund.age = NULL,catch.age = NULL,preload=T)
{
  
  if(preload==F)
  {
  # Download the functions we'll need from github
  funs <- c("https://raw.githubusercontent.com/dave-keith/ICM/main/Scripts/functions/simple_Lotka_r.r",
            "https://raw.githubusercontent.com/dave-keith/ICM/main/Scripts/functions/backwards_project.r",
            "https://raw.githubusercontent.com/dave-keith/ICM/main/Scripts/functions/forward_project.r"
            
  )
  # Now run through a quick loop to load each one, just be sure that your working directory is read/write!
  for(fun in funs) 
  {
    download.file(fun,destfile = basename(fun))
    source(paste0(getwd(),"/",basename(fun)))
    file.remove(paste0(getwd(),"/",basename(fun)))
  }
  } # end if(preload==F)
  
  # Uncomment these lines if troubleshooting, you'll need to change the directories.
  #source("D:/Github/ICM/Scripts/functions/simple_Lotka_r.r")
  #source("D:/Github/ICM/Scripts/functions/forward_project.r")
  #source("D:/Github/ICM/Scripts/functions/backwards_project.r")
  
  require(optimx)  || stop("Please load the 'optimx' package which you'll need for the optimations to run")
  
 
  #Initialize a bunch of objects
  n.years<-length(years)
  # Not using this, but I might someday, this is the minimum age in the VPA, we've filled in the rest of the 0s...
  min.age <- min(ages)
  #  Run from min age to max age!
  ages <- min.age:max(ages)
  # Making a dataframe for later.
  res <- data.frame(year = years,vpa.abund = N.init, est.abund = NA, diff= NA,per.diff = NA,lambda=NA,lambda.vpa.init = NA,
                    removals = NA,removals.init = NA,mean.fec=NA,mean.vpa.fec = NA,mean.nm = NA,mean.vpa.nm = NA)
  # making sure objects are nicely formated
  fecund.ts <- data.frame(fecund)
  nm.ts <- data.frame(nm)
  fm.ts <- data.frame(fm)
  z.ts <- data.frame(fm.ts + nm.ts)
  catch.age.ts <- data.frame(catch.age)
  abund.age.ts <- data.frame(abund.age)
  
  # Get the average removals so we can identify really low years
  tot.catch <- rowSums(catch.age.ts,na.rm=T)
  mn.fm.age <- colMeans(fm.ts,na.rm=T)
  mn.catch <- mean(tot.catch,na.rm=T)
  
  # an index to move through the years
  if(direction == 'backwards') y.index <- (n.years-1):1 # Starting at the second last year of data makes the indexing below way more consistent between directions
  if(direction == 'forwards')  y.index <- 1:(n.years-1) # Stop at the second last year of data
  # Get the initial abundance
  N.up <- N.init
  # Now move through the years making the model fit the VPA results.
  for(y in y.index)
  {
    print(paste0("Year = ",years[y]))
    pop.est <- N.up
    # The first step is to fit the default lotka.r values
    # When going forwards we use everything in the 'current' year to project to next year
    # When going backwards we use everything in the 'previous' year to project from the next year
    # Getting the data we need to start the year
    nm.tmp <- nm.init <- nm[y,]
    fecund.tmp <- fecund.init <-  fecund[y,]
    fm.tmp <- fm.init <- fm[y,]
    z.tmp <- z.init <- nm.tmp + fm.tmp
    catch.age.tmp <- catch.age.init <- catch.age[y,]
    tot.catch.tmp <- tot.catch.init <- sum(catch.age.tmp)
    abund.age.tmp <- abund.age.ts[y,]
    
    # If claiming no fishing, make there be a little bit of fishing so we have something to tweak if we want to.
    # I've set it to be 1% of the mean fishing mortality, so this will be tiny (F of 0.4 would be 0.004)
    # FIX: Do we need to do this?
    if(max(fm.tmp,na.rm=T) == 0)
    {
      fm.tmp <- 0.01*mn.fm.age 
      catch.age.tmp <- (1-exp(-fm.tmp)) * abund.age.tmp
      tot.catch.tmp <- sum(catch.age.tmp,na.rm=T)
    }
    
    # FIX: So for our removals time series, we put the removals between y+1 and y 
    # down as year y+1.  We can change this, but that's how this is set up at the moment.
    # So when going forwards we want to use the current year index, but when going backwards
    # we want to use the previous years. This is all handled by the y.index happily.
    # Going to line up the r with the 'previous' year, so the r in column 2020 is the r used to move between 2020 and 2021 (in either direction)
    
    # Here's is where the magic happens and we get our estiamte of population growth for the year based on the mortality and fecundity info
    lotka.est <- simple.lotka.r(mort = z.tmp,fecund=fecund.tmp,ages=ages)
    # Pull out the estimate from this
    #browser()
    r.est = lotka.est$res$r
    # If the r.est is bonkers, we reduce fecundity until it isn't
    while(r.est > 10) 
    {
      # if the lotka isn't converging set the age class fecundity to be 90% of what we had if tuning on fecundity
      if(tuner %in% c("fecund","fec_nm")) fecund.tmp <- 0.9*fecund.tmp 
      # if tuning on natural mortality increase that by 10%
      if(tuner %in% c("nm","fec_nm",'z')) nm.tmp <- nm.tmp*1.1
      # if tuning on z, increase fishing mortality by 10% too.
      if(tuner == "z") fm.tmp <- fm.tmp*1.1
      # If chaning nm or fm then we need a new z.tmp object that accounts for that
      if(tuner != 'fecund') z.tmp <- nm.tmp + fm.tmp
      
      lotka.est <- simple.lotka.r(mort = z.tmp,fecund=fecund.tmp,ages=ages)
      # Pull out the estimate from this
      #browser()
      r.est = lotka.est$res$r
      print("lotka wasn't converging, we pre-tuned the tuner variables in steps of 10% until lotka converged")
    }
 
    # If going backwards we do the back.projection
    if(direction == 'backwards')
    {
        # Do the projection, using the exponential model (because we have the 'real' estimate of r)
        exp.res <- back.proj(option = "exponential",pop.next = pop.est,r=r.est)
        # get the abundance
        pop.step <- exp.res$Pop.current
        # calculate the difference and % difference
        diff.org <- abund.ts[y] - pop.step
        per.diff.org <- 100*((pop.step - abund.ts[y])/abund.ts[y])
    } # end if(direction == 'backwards')
    
    # If going forwards we do the forward projection
    if(direction == 'forwards')
    {
      #browser()
        # Do the projection, using the exponential model (because we have the 'real' estimate of r)
        exp.res <- for.proj(option = "exponential",pop.last = pop.est ,r=r.est)
        if(exp.res$Pop.current < 0) exp.res$Pop.current =0 # don't let it drop below 0
        # Update the abundance
        pop.step <- exp.res$Pop.current
        # calculate the difference and % difference
        diff.org <- abund.ts[y+1] - pop.step
        per.diff.org <- 100*((pop.step - abund.ts[y+1])/abund.ts[y+1])
    } # end if(direction == 'forwards')
    # rename the diff and % difference objects
    diff <- diff.org
    per.diff <- per.diff.org
    
    
    # So now we compare that value to the 'real' value from the vpa
    #Initial conditions
    count =0
    last.per.diff <- 1e6 # Make it ridiculous high to make sure we get started...
    # Run this until one of these becomes false.
    get.better <- T    # get.better is making sure we are going the right way
    big.diff <- T      # big.diff means we are far away so we can make a big change
    not.close.enough <- T  # Not close enough means we are too far away from correct so we keep going.
    # So as long as these three things are TRUE we are in this while. This while is the 'it's a big difference' while
    while(get.better && big.diff &&  not.close.enough) 
    {
      #if(y %in% 35:37) browser()
      if(last.per.diff == per.diff & count > 2) stop("You broke something the percent differences aren't changing.")
      # Update
      last.per.diff <- per.diff
      # indexing...
      count = count + 1
      # Depending on how you decided to 'tune' depends which if we pop into
      # For the first step we need to figure out what direction to go, thus for count = 1 we check decline and increase and have a decently big step size.
      if(count == 1)
      {
        # if varying fecundity do this...
        if(tuner %in% c('fecund',"fec_nm"))
        {
          fecund.inc <- fecund.tmp * (1+5*step.size)
          fecund.dec <- fecund.tmp / (1+5*step.size)
        } # end if tuner is fecund or both
        
        # If varying natural mortality go in here, note that 
        if(tuner %in% c('nm','z',"fec_nm"))
        {
          nm.inc <- nm.tmp * (1+5*step.size)
          nm.dec <- nm.tmp / (1+5*step.size)
          if(tuner %in% c('nm',"fec_nm"))
          {
            z.inc <- nm.inc + fm.tmp
            z.dec <- nm.dec + fm.tmp
          } # end if tuner is nm or both
          } # end if tuner is nm, z, or both
          
          # If varying f go in here
          if(tuner %in% c('fm','z'))
          {
            fm.dec <- fm.tmp / (1+5*step.size)
            catch.age.dec <- (1-exp(-fm.dec)) * abund.age.tmp
            catch.dec <- sum(catch.age.tmp,na.rm=T)
            fm.inc <- fm.tmp * (1+5*step.size)
            catch.age.inc <- (1-exp(-fm.inc)) * abund.age.tmp
            catch.inc <- sum(catch.age.tmp,na.rm=T)
            # if only tweaking fishing mortality
            if(tuner == 'fm')
            {
            z.inc <- nm.tmp + fm.inc
            z.dec <- nm.tmp + fm.dec 
            } # end the fishing mortality tuner 
          } #end tuner if varying fishing mortality
          
          # If tweaking both f and natural mortality (note this grabs the natural mortality and fishing mortalities from the above if loops)
          if(tuner == 'z')
          {
            z.inc <- nm.inc + fm.inc
            z.dec <- nm.dec + fm.dec 
          } # end tuner is total mortality

          # Because we are in count 1, we calculate the lotka both directions and figure out which way is making the fit better...
          # what we run depends on our tuner...
          if(tuner %in% c('fm','nm','z'))
          {
            lotka.inc <- simple.lotka.r(mort = z.dec,fecund=fecund.tmp,ages=ages)
            lotka.dec <- simple.lotka.r(mort = z.inc,fecund=fecund.tmp,ages=ages)
          } # end fm/nm tuner
            
          if(tuner == 'fecund') 
          {
            lotka.inc <- simple.lotka.r(mort = z.tmp,fecund=fecund.inc,ages=ages)
            lotka.dec <- simple.lotka.r(mort = z.tmp,fecund=fecund.dec,ages=ages)
          } # end fecund tuner
            
          if(tuner == 'fec_nm') 
          {
            lotka.inc <- simple.lotka.r(mort = z.dec,fecund=fecund.inc,ages=ages)
            lotka.dec <- simple.lotka.r(mort = z.inc,fecund=fecund.dec,ages=ages)
          } # end fecund & nm tuner  
            
          r.inc <- lotka.inc$res$r
          r.dec <- lotka.dec$res$r
          
          # If doing a backwards simulation we see how close we got to the 'vpa' abundance using the back project function
          if(direction == 'backwards')
          {
            # Get your population estimates by increasing r and decreasing r
            exp.inc <- back.proj(option = "exponential",pop.next = pop.est,r=r.inc)
            pop.est.inc <- exp.inc$Pop.current
            exp.dec <- back.proj(option = "exponential",pop.next = pop.est,r=r.dec)
            pop.est.dec <- exp.dec$Pop.current
            # Now we test the two methods...
            diff.inc <- abund.ts[y] - pop.est.inc
            per.diff.inc <- 100*((pop.est.inc - abund.ts[y])/abund.ts[y])
            diff.dec <- abund.ts[y] - pop.est.dec
            per.diff.dec <- 100*((pop.est.dec - abund.ts[y])/abund.ts[y])
          } # end if(direction == 'backwards')
          
          # If doing a forwards simulation we see how close we got to the 'vpa' abundance using the forwards project function
          if(direction == 'forwards')
          {
            # Get your population estimates by increasing r and decreasing r
            exp.inc <- for.proj(option = "exponential",pop.last = pop.est ,r=r.inc)
            if(exp.inc$Pop.current < 0) exp.inc$Pop.current =0 # don't let it drop below 0
            pop.est.inc <- exp.inc$Pop.current
            # And the decrease one
            exp.dec <- for.proj(option = "exponential",pop.last = pop.est ,r=r.dec)
            if(exp.dec$Pop.current < 0) exp.dec$Pop.current =0 # don't let it drop below 0
            pop.est.dec <- exp.dec$Pop.current
            # Now we test the two methods...
            diff.inc <- abund.ts[y+1] - pop.est.inc
            per.diff.inc <- 100*((pop.est.inc - abund.ts[y+1])/abund.ts[y+1])
            diff.dec <- abund.ts[y+1] - pop.est.dec
            per.diff.dec <- 100*((pop.est.dec - abund.ts[y+1])/abund.ts[y+1])
          } # end if(direction == 'forwards')
        
       
        # we are getting better by increasing r, i.e. either increasing fecundity or decreasing F or M
        if(abs(per.diff.inc) < abs(per.diff.dec))
        {
          going <- 'up' # up means increasing r
          # now pick the 'correct' new values of the tuned parameters as appropriate
          if(tuner %in% c('fecund',"fec_nm")) fecund.tmp <- fecund.inc
          if(tuner %in% c('nm','z',"fec_nm"))     nm.tmp <- nm.dec
          if(tuner %in% c('fm','z'))    
          {
            fm.tmp <- fm.dec
            catch.age.tmp <- catch.age.dec
          }  
          per.diff <- per.diff.inc

        } else 
          {
            going <- 'down' # down means decreasing r
            # now pick the 'correct' new values of the tuned parameters as appropriate
            if(tuner %in% c('fecund',"fec_nm")) fecund.tmp <- fecund.dec
            if(tuner %in% c('nm','z',"fec_nm"))    nm.tmp <- nm.inc
            if(tuner %in% c('fm','z'))    
            {
              fm.tmp <- fm.inc
              catch.age.tmp <- catch.age.inc
            }
            per.diff <- per.diff.dec
            #pop.tunes <- pop.est.dec
          } # end if else
        } # end count ==1
        
        # OK, so we have our direction figured out, now we go off in that direction and optimize
        if(count > 1)
        {
          # Increasing the population growth rate
          if(going == 'up')  
          {
            if(tuner %in% c('fecund',"fec_nm")) fecund.tmp <- fecund.tmp*(1+step.size)
            if(tuner %in% c('nm','z',"fec_nm")) nm.tmp <- nm.tmp / (1+step.size)
            if(tuner %in% c('fm','z')) 
            {
              fm.tmp <- fm.tmp / (1+step.size)
              catch.age.tmp <- (1-exp(-fm.tmp)) * abund.age.tmp
              tot.catch.tmp <- sum(catch.age.tmp,na.rm=T)
            } # end if(tuner == 'fm') 

            # So once we have the ball rolling, let's see how close to correct we are
            if(count > 2)
            {
              # When the percent difference is > 5% we can take some relatively big steps to get closer
              if(abs(per.diff) > 5) 
              {
                # Go for big steps and for everybody, the numbers here seem arbitrary, but I set them to values that 
                # quickly resulted in convergence without overshooting the optimum.
                fec.ss <- nm.ss <- fm.ss <- 5*step.size
                if(abs(per.diff)>15) fec.ss <- nm.ss <- fm.ss <- 3*fec.ss
                if(tuner %in% c('fecund',"fec_nm")) fec.ss <- 10*step.size
                if(tuner %in% c('nm','z',"fec_nm")) nm.ss <- 10*step.size
                if(tuner %in% c('fm','z'))  fm.ss <- 10*step.size
                
                fecund.tmp <- fecund.tmp*(1+(fec.ss))  # if(tuner =='fecund')
                nm.tmp <- nm.tmp/(1+(nm.ss))
                if(tuner != c("fec_nm"))
                {
                  fm.tmp <- fm.tmp / (1+(fm.ss))
                  catch.age.tmp <- (1-exp(-fm.tmp)) * abund.age.tmp
                  tot.catch.tmp <- sum(catch.age.tmp,na.rm=T)
                } # End don't do this if we just want to tune on fecundity and natural mortality
              }# end abs(per.diff > 5)
            } # end if count > 2
          } # end if going is up.
          
          
          # Now if we need to decrease the population growth rate
          if(going == 'down') 
          {
            if(tuner %in% c('fecund',"fec_nm")) fecund.tmp <- fecund.tmp/(1+step.size)
            if(tuner %in% c('nm','z',"fec_nm")) nm.tmp <- nm.tmp*(1+step.size)
            if(tuner %in% c('fm','z'))
            {
              fm.tmp <- fm.tmp * (1+step.size)
              catch.age.tmp <- (1-exp(-fm.tmp)) * abund.age.tmp
              tot.catch.tmp <- sum(catch.age.tmp,na.rm=T)
            } # end if(tuner == 'fm')
            # Once r is big, we also turn down nm!!
            if(count > 2)
            {
              # Go for big steps and for everybody, the numbers here seem arbitrary, but I set them to values that 
              # quickly resulted in convergence without overshooting the optimum.
              if(abs(per.diff) > 10) 
              {
                fec.ss <- nm.ss <- fm.ss <- 2*step.size
                if(tuner %in% c('fecund',"fec_nm")) fec.ss <- 10*step.size
                if(tuner %in% c('nm','z',"fec_nm")) nm.ss <- 10*step.size
                if(tuner %in% c('fm','z')) 
                {
                  fm.ss <- 10*step.size
                  # Really far away means really big jump, if fishing mortality is super low we bump it up a bunch
                  # only if we are tuning on z and fishing mortality
                  if(abs(per.diff) > 15)
                  {
                    fm.big.jump <- tot.catch.tmp < 0.1*mn.catch
                    if(fm.big.jump) fm.tmp <- 0.5*mn.fm.age else fm.ss <- 20*step.size
                  } # end the bigger than 15% jump...
                } # end the fm and z 
                # Decrease fecundity when we need a big decrease in R 
                fecund.tmp <- fecund.tmp/(1+fec.ss)  # if(tuner =='fecund')
                nm.tmp <- nm.tmp*(1+(nm.ss))
                if(tuner != c("fec_nm"))
                {
                  fm.tmp <- fm.tmp * (1+fm.ss)
                  catch.age.tmp <- (1-exp(-fm.tmp)) * abund.age.tmp
                  tot.catch.tmp <- sum(catch.age.tmp,na.rm=T)
                } # end the don't do it if only tuning on fec and nm.
              }# end abs(per.diff > 30)
            } # end if count > 2
          } # end going down 
          
          # the the total Z from the above, use that and fecundity and get the new lotka r estimate.
          z.tmp <- nm.tmp + fm.tmp
          lotka.tunes <- simple.lotka.r(mort = z.tmp,fecund=fecund.tmp,ages=ages)
          r.tunes <- lotka.tunes$res$r
        
        # Run this throught the backwards or forwards simulations as appropriate
        if(direction == 'backwards')
        {
          exp.tunes <- back.proj(option = "exponential",pop.next = pop.est,r=r.tunes)
          pop.tunes <- exp.tunes$Pop.current
          diff.tunes <- abund.ts[y] - pop.tunes
          per.diff.tunes <- 100*((pop.tunes - abund.ts[y])/abund.ts[y])
          per.diff <- per.diff.tunes
        } # end if(direction == 'backwards')
        
        # And forwards
        if(direction == 'forwards')
        {
          exp.tunes <- for.proj(option = "exponential",pop.last = pop.est ,r=r.tunes)
          if(exp.tunes$Pop.current < 0) exp.tunes$Pop.current =0 # don't let it drop below 0
          pop.tunes <- exp.tunes$Pop.current
          diff.tunes <- abund.ts[y+1] - pop.tunes
          per.diff.tunes <- 100*((pop.tunes - abund.ts[y+1])/abund.ts[y+1])
          per.diff <- per.diff.tunes
        } # end if(direction == 'forwards')
     
      } # end if count >1
        print(paste0("Percent difference = ",round(per.diff))) 
        
       
        # Now we look at the results of the above to see if we can exist the while 'loop'.
        # Basically we are looking to see if the simulations are improving, 
        # if we are 'stuck' and if we are close enough to right to bail out of the while
        # Get out of the while if the results are getting worse (we've gone past the optimum result)
        if(abs(per.diff) > abs(last.per.diff) ) get.better <- F
        # Get out if the difference is less than 1% of the step size and we have run 3 steps, basically
        # avoids infinite loops (stop if it doesn't do anything anymore)
        if((abs(abs(per.diff) - abs(last.per.diff))) < step.size/100 & count > 2) big.diff <- F
        # If the % difference is < 5%, the absolute percent difference is less than 20 times our step size, and we have run through 3 steps
        # Basically we are 'close enough
        # FIX: 5% is arbitrary, we can change this, probably makes sense to be a function option.
        if(abs(per.diff)  < 5 & abs(per.diff - last.per.diff) < step.size*20 & count > 2) not.close.enough <- F
      } # end while(get.better && big.diff &&  not.close.enough)
      
      # If the maximum count = 1 then the original is the best fit!
      if(max(count) == 1) 
      {
          r.tunes <- r.est
          pop.tunes <- pop.step
          diff.tunes <- diff.org
          per.diff <- per.diff.org
          fecund.tmp <- fecund.init
          nm.tmp <- nm.init
      } # end if(max(count) == 1) 
      
    # Get the results tidied up when going forward, this is all about indices.
    if(direction == 'forwards')  
    {
      index <- which(years == years[y+1])
      #browser() 
      res$lambda[index-1] <- exp(r.tunes)
      #browser()
      res$lambda.vpa.init[index-1] <- exp(r.est)
      res$removals[index-1] <- tot.catch.tmp
      res$removals.init[index-1] <- tot.catch.init
      res$mean.fec[index-1] <- mean(as.numeric(fecund.tmp))
      res$mean.vpa.fec[index-1] <- mean(as.numeric(fecund.init))
      res$mean.nm[index-1] <- mean(as.numeric(nm.tmp))
      res$mean.vpa.nm[index-1] <- mean(as.numeric(nm.init))
      fecund.ts[index-1,] <- fecund.tmp
      nm.ts[index-1,] <- nm.tmp
      fm.ts[index-1,] <- fm.tmp
      z.ts[index-1,] <- z.tmp
    } # end forwards
    # FIX... These indicies might be wrong for going backwards now
    if(direction == 'backwards')
    {
      index <- which(years == years[y])
      res$lambda[index+1] <- exp(r.tunes)
      res$lambda.vpa.init[index+1] <- exp(r.est)
      res$removals[index+1] <- tot.catch.tmp
      res$removals.init[index+1] <- tot.catch.init
      res$mean.fec[index+1] <- mean(as.numeric(fecund.tmp))
      res$mean.vpa.fec[index] <- mean(as.numeric(fecund.init))
      res$mean.nm[index+1] <- mean(as.numeric(nm.tmp))
      res$mean.vpa.nm[index] <- mean(as.numeric(nm.init))
      fecund.ts[index+1,] <- fecund.tmp
      nm.ts[index+1,] <- nm.tmp
    }
 
    # Summary stuff
    res$vpa.abund[index] <- abund.ts[index]
    res$est.abund[index] <- pop.tunes
    res$diff[index] <- diff.tunes
    res$per.diff[index] <- per.diff
 
    # Get the next N ready...
    N.up <- pop.tunes
    
    # We then save that result and re-run everything with a new value of fecundity or natural mortality, if that value improves our
    # estimate we keep doing that.  Once it no longer improves the value we stop, save the value, and use that population value to start 
    # the next round.
    # I think what we can have the 'both' scenario do is to change fecundity in the first iteration
    # then the second iteration we change mortality, the third iteration changes fecundity... etc etc.
    # Then we need to break out of the iterative loop when the result gets worse (in terms of absolute % difference)
    # Save that result then move to the next year.... 
    #browser() 
  } # end the years loop
  
  return(list(res=res,fecund.opt = fecund.ts,nm.opt = nm.ts,fm.opt = fm.ts,z.opt = z.ts,fecund.org = fecund,nm.org=nm,fm.org = fm,z.org = nm+fm))
} # end the function  
  
  
  
  
  
  
  
  
  
  