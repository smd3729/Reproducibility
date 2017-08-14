remove(list=ls())
library(skellam)
library(xlsx)
library(bbmle)

resultFilename = 'james_record_v9.csv'
ciFilename     = 'james_record_v9_CI.csv'


########################################################
# Get the data
#
########################################################
data         = read.xlsx('Reproducibility dataset 2.xlsx', sheetIndex=2)


maxitLoc     = 50000
methodLoc    = 'Nelder-Mead'
optimizerLoc = 'optim'

# Basic parameter list
# beta0 = intercept
# beta1 = NumPubs, continuous
# beta2 = pubs/year, continuous
# beta3 = Yrs1stpub, continuous
# beta4 = Country, catagorical N levels
# beta4a = Region, Latin America or not Latin America
# beta5 = Publspeciesdescript, catagorical N levels
# beta5a = Expert, generists vs specialist 
# beta6 = Species concept, catagorical N levels
# beta7 = Taxon Assigned, catagorical N levels
# beta7a = #spp/genus, continuous
# beta7b = %genus, continuous
# beta7c = Genus age, continuous
# beta7d = GShits, continuous
# beta7e = File KB, continuous
# beta7f = Speciems, continuous
# beta8 = DataTypeUsed, catagorical N levels
# beta9 = MorphSampleSufficient, catagorical N levels
# beta10 = MolecSampleSufficient, catagorical N levels
# beta11 = Useful, catagorical N levels
# beta12 = ConsultLit, catagorical N levels
# beta13 = LitInfluence, catagorical N levels
# beta14 = NumTasks, ordinal
# beta14a = Allignment, catagorical 2 levels
# beta14b = Geneticdistance, catagorical 2 levels
# beta14c = Modeltesting, catagorical 2 levels
# beta14d = Multivariate stats, catagorical 2 levels
# beta14e = Phylogenetics, catagorical 2 levels
# beta14f = Simplestats, catagorical 2 levels
# beta14g = AutoSpeciesdelim, catagorical 2 levels
# beta15  = AutomDiff

# beta16 and higher are interaction terms. 

########################################################
# Define the model and negLogL function
#
########################################################
model = function(params, data)
{
    yHat =  params['beta0'] + # Intercept

            params['beta1']     * data$NumPubs                 + #NumPubs

            #params['beta2']     * data$pubs.year                 + # pubs/year

            params['beta4_Col'] * ( data$Country == 'Colombia')  + # Country
            params['beta4_Mex'] * ( data$Country == 'Mexico')    + # Country
            params['beta4_Swe'] * ( data$Country == 'Sweden')    + # Country
            params['beta4_USA'] * ( data$Country == 'USA')       + # Country
            params['beta4_Ven'] * ( data$Country == 'Venezuela') + # Country

            #params['beta4a_NotLA'] * ( data$Region == 'NotLA' ) +

            params['beta5_Aca'] * (data$Publspeciesdescrip == 'Acanthomorpha' ) +  #Publspeciesdescript
            params['beta5_All'] * (data$Publspeciesdescrip == 'All'           ) +  #Publspeciesdescript
            params['beta5_Cyp'] * (data$Publspeciesdescrip == 'Cypriniformes' ) +  #Publspeciesdescript
            params['beta5_Gym'] * (data$Publspeciesdescrip == 'Gymnotiformes' ) +  #Publspeciesdescript
            params['beta5_Sil'] * (data$Publspeciesdescrip == 'Siluriformes'  ) +  #Publspeciesdescript

            #params['beta5a_Specialist'] * (data$Expert == 'Specialist'  ) +  #Expert

            params['beta6_Mor'] * (data$Species.concept == 'Morphophenetic')    +  #Species concept
            params['beta6_Nom'] * (data$Species.concept == 'Nominalist')        +  #Species concept
            params['beta6_Phy'] * (data$Species.concept == 'Phylogenetic')      +  #Species concept
            params['beta6_Rep'] * (data$Species.concept == 'Reproductive')      +  #Species concept

            params['beta7_Neo'] * ( data$Taxon.Assigned == 'Neoplecostomus')    +  #Taxon Assigned
            params['beta7_Pro'] * ( data$Taxon.Assigned == 'Profundulus')       +  #Taxon Assigned
            params['beta7_Tet'] * ( data$Taxon.Assigned == 'Tetragonopterus')   +  #Taxon Assigned

            #params['beta8_MOR'] * ( data$DataTypeUsed == 'MORPH' )              +  #DataTypeUsed
            #params['beta8_MOL'] * ( data$DataTypeUsed == 'MOLEC' )              +  #DataTypeUsed

            params['beta11_MOR'] * ( data$Useful == 'MORPH' )              +  #DataTypeUsed
            params['beta11_MOL'] * ( data$Useful == 'MOLEC' )              +  #DataTypeUsed

            params['beta14'] * data$NumTasks + #NumTasks

            #params['beta14a'] * (data$Allignment         == 'YES' ) +
            #params['beta14b'] * (data$Geneticdistance    == 'YES' ) +
            #params['beta14c'] * (data$Modeltesting       == 'YES' ) +
            #params['beta14d'] * (data$Multivariate.stats == 'YES' ) +
            #params['beta14e'] * (data$Phylogenetics      == 'YES' ) +
            #params['beta14f'] * (data$Simplestats        == 'YES' ) +
            #params['beta14g'] * (data$AutoSpeciesdelim   == 'YES' ) 

            params['beta15'] * data$AutomDiff #Difference in scientists score between automation of species deliniation

    return(yHat)
}

dskellam3 = function(x, mu, shape, log=F)
{
    sigma2 = max(abs(mu)) + shape
    return( dskellam(x, lambda1 = (sigma2 + mu)/2.0, lambda2 = (sigma2 - mu)/2.0, log = log) )
}

shape_map = function( x ) 
{
    if ( x < 0  ) { return( exp(x) ) }
    return( x + 1 )
}
shape_map_inv = function( y ) 
{
    if ( y < 1 ) { return ( log(y) ) }
    return ( y - 1 )
}

negLogL_for_mle2 = function( params, data )
{
    paramsLoc          = params
    #paramsLoc['shape'] = exp(params['shape'])
    paramsLoc['shape'] = shape_map(params['shape'])
    yHat               = model( paramsLoc, data)
    ret                = -sum( dskellam3( data$Deviation, mu = yHat, shape=paramsLoc['shape'], log=T) )
    return(ret)
}

parnames(negLogL_for_mle2) = c( 'beta0',
                               'beta1',
                               'beta4_Col',
                               'beta4_Mex',
                               'beta4_Swe',
                               'beta4_USA',
                               'beta4_Ven',
                               'beta5_Aca',
                               'beta5_All',
                               'beta5_Cyp',
                               'beta5_Gym',
                               'beta5_Sil',
                               'beta6_Mor',
                               'beta6_Nom',
                               'beta6_Phy',
                               'beta6_Rep',
                               'beta7_Neo',
                               'beta7_Pro',
                               'beta7_Tet',
                               'beta11_MOR',
                               'beta11_MOL',
                               'beta14',
                               'beta15',
                               'shape' )


########################################################
# Define some analysis functions
#
########################################################

profile2 = function( mle2Fit, varName )
{
    cat('profile2: Msg: Extracting information\n')
    dataLoc    = mle2Fit@data
    varEst     = mle2Fit@fullcoef[varName]
    paramFull  = mle2Fit@fullcoef
    nllFull    = negLogL_for_mle2( paramEst, data)

    ret             = data.frame ( matrix( mle2Fit@fullcoef, nrow=1) )
    names(ret)      = names( mle2Fit@fullcoef )
    ret$convergence = mle2Fit@details$convergence
    ret$negLogL     = mle2Fit@min
    ret$delta       = ret@negLogL - nullFull
    ret$chisq       = 2*( ret$delta )
    ret$p           = pchisq( ret$chisq, df = 1 )
    ret$type        = 'Fit'

    cat('profile2: Msg: Step 1 computing slice\n')
    # Step 1: Compute a slice to get an estimate on what
    # range of values to use.
    delta = abs( 0.95 * varEst )
    sliceMin = varEst
    p = 0
    while ( p < 0.99 ) {
        sliceMin          = sliceMin - delta

        paramEst          = paramFull
        paramEst[varName] = sliceMin
        nllSlice          = negLogL_for_mle2( paramEst, data)
        chiSqSlice        = 2*( nllSlice - nllFull )
        p                 = pchisq(chiSqSlice, df=1)
        #cat('profile2: Msg:\n')
        #cat('\tsliceMin = ',sliceMin,'\n')
        #cat('\tnllSlice = ',nllSlice,'\n')
        #cat('\tnllFull  = ',nllFull,'\n')
        #cat('\tchiSqEst = ',chiSqEst,'\n')
        #cat('\tp        = ',p,'\n')
    }

    profileMin = sliceMin 
    p          = 0
    while ( p < 0.99 ) {
        paramEst          = paramFull
        profileMin        = profileMin - delta
        paramEst[varName] = profileMin
        paramFixed        = c( profileMin )
        names(paramFixed) = varName

        fitFixed          = mle2( minuslogl = negLogL_for_mle2, 
                                  start     = paramEst, 
                                  method    = methodLoc,
                                  optimizer = optimizerLoc,
                                  fixed     = paramFixed,
                                  data      = list( data = data),
                                  control   = list( maxit = maxitLoc) )

        nllFixed          = netgLogL_for_mle2( fitFixed$fullcoef, data ) 
        chiSqFixed        = 2*( nllFixed - nllFull )
        p                 = pchisq( chiSqFixed, df = 1)
    }

    sliceMax = varEst
    p = 0
    while ( p < 0.99 ) {
        sliceMax          = sliceMax + delta
        paramEst[varName] = sliceMax
        nllSlice          = negLogL_for_mle2( paramEst, data)
        chiSqEst          = 2*( nllSlice - nllFull)
        p                 = pchisq(chiSqEst, df=1)
        #cat('profile2: Msg:\n')
        #cat('\tsliceMin = ',sliceMin,'\n')
        #cat('\tnllSlice = ',nllSlice,'\n')
        #cat('\tnllFull  = ',nllFull,'\n')
        #cat('\tchiSqEst = ',chiSqEst,'\n')
        #cat('\tp        = ',p,'\n')
    }

    profileMax = sliceMax
    p          = 0
    while ( p < 0.99 ) {
        profileMax        = profileMax + delta
        paramEst[varName] = profileMax
        paramFixed        = c( profileMax )
        names(paramFixed) = varName

        fitFixed          = mle2( minuslogl = negLogL_for_mle2, 
                                  start     = paramEst, 
                                  method    = methodLoc,
                                  optimizer = optimizerLoc,
                                  fixed     = paramFixed,
                                  data      = list( data = data),
                                  control   = list( maxit = maxitLoc) )

        nllFixed          = netgLogL_for_mle2( fitFixed$fullcoef, data ) 
        chiSqFixed        = 2*( nllFixed - nllFull )
        p                 = pchisq( chiSqFixed, df = 1)
    }


    cat('profile2: Msg: Step 2 Creating profile range\n')
    # Step 2: Create a range of values to test
    #ret          = data.frame( c( varValue, seq(sliceMin, sliceMax, by=delta * 0.1)) )
    profile          = c( varEst, seq(profileMin, profileMax, length.out=14))
    profile          = profile[ order(profile) ]

    cat('profile2: Msg: Step 3 Computing profile\n')
    # Step 3: Compute the profile
    for ( i in 1:length(profile) ) {
        cat('profile2: Msg: i = ',i,' of ',nrow(profile),'\n')
        value              = profile[i]

        paramEst           = paramFull
        paramEst[varName]  = value
        fixedLoc           = c( value )
        names(fixedLoc)    = varName

        fitLoc             = mle2( minuslogl = negLogL_for_mle2, 
                                   start     = paramEst, 
                                   method    = methodLoc,
                                   optimizer = optimizerLoc,
                                   fixed     = fixedLoc,
                                   data      = list( data = data),
                                   control   = list(maxit=maxitLoc) )

        nllProfile          = negLogL_for_mle2( fitLoc@fullcoef, data)
        delta               = nllProfile - nllFull
        chiSqEst            = 2*( delta )
        p                   = pchisq( chisqEst, df = 1)

        newRow              = data.frame( matrix( fitLoc@fullcoef, nrow=1 ) )
        names(newRow)       = names(fitLoc@fullcoef)
        newRow$convergence  = fitLoc@details$convergence
        newRow$negLogL      = fitLoc@min
        newRow$delta        = delta
        newRow$chisq        = chiSqEst
        newRow$p            = p
        newRow$type         = varName

        ret                = rbind( ret, newRow )
    }

    cat('profile2: Msg: Step 5: Done\n')
    return(ret)
}

confint2 = function( result ) 
{
    paramList = names(result@fullcoef)
    paramList = paramList[ !paramList%in%c('shape') ]
    ret       = data.frame()
    for ( i in 1:length(paramList) ) {
        parName   = paramList[i]
        resultCI2 = confint2_var( result, parName )

        cat(parName, resultCI2$lowerCI[1], resultCI2$upperCI[1], sep=',')
        cat('\n')

        if ( nrow(ret) == 0 ) {
            ret = resultCI2
        } else {
            ret = rbind( ret, resultCI2 )
        }
    }
    return(ret)
}

confint2_var = function( fit, varName, alphaLevel=0.95, threshold=0.01, stepScale=0.1, maxIt = 50, verbose = F )
{
    # Step 1: Gather the baisc information
    data     = fit@data
    if ( class(data) == 'list' && length(data) == 1 && class(data[[1]]) == 'data.frame' ) {
        data = data[[1]]
    }
    if ( class(data) != 'data.frame' ) {
        cat('pois.R: profile3: Fuck R\n')
        return(0)
    }
    paramEst = fit@fullcoef
    nllFull  = fit@min

    # Step 2: Do a modified binary search for the lower CI
    varEst       = paramEst[ varName ]
    stepSize     = stepScale * abs( varEst )
    delta        = - stepSize
    lowerCI      = varEst
    p            = 0
    count        = 0
    lowerP       = 0
    lowerChiSq   = 0
    lowerConverg = NA

    while ( abs( p - alphaLevel ) > threshold && count < maxIt ) {
        count             = count + 1

        paramGuess        = paramEst

        lowerCI           = lowerCI + delta
        paramFixed        = c( lowerCI )
        names(paramFixed) = varName

        fitFixed = mle2( minuslogl = negLogL_for_mle2,
                        start = paramGuess,
                        fixed = paramFixed,
                        data  = list( data = data) )

        nllFixed = fitFixed@min
        chiSq    = 2*( nllFixed - nllFull)
        p        = pchisq( chiSq, df = 1)

        lowerChiSq   = chiSq
        lowerP       = p
        lowerConverg = fitFixed@details$convergence

        if ( verbose ) {
            cat('lowerCI = ',lowerCI,'\n')
            cat('chiSq   = ',chiSq,  '\n')
            cat('p       = ',p,      '\n')
            cat('delta   = ',delta,  '\n\n')
        }


        # if ( p < alphaLevel && delta < 0 ) { } # keep going

        if ( p > alphaLevel && delta < 0 ) { delta = - delta/2.0 } # switch directions

        #if ( p > alphaLevel && delta > 0 ) { } # keep going 

        if ( p < alphaLevel && delta > 0 ) { delta = - delta/2.0 } #switch directions

    }

    


    # Step 3: Do a modified binary search for the lower CI
    varEst       = paramEst[ varName ]
    stepSize     = stepScale * abs( varEst )
    delta        = stepSize
    upperCI      = varEst
    p            = 0
    count        = 0
    upperP       = 0
    upperChiSq   = 0
    upperConverg = NA

    while ( abs( p - alphaLevel ) > threshold && count < maxIt ) {
        count             = count + 1

        paramGuess        = paramEst

        upperCI           = upperCI + delta
        paramFixed        = c( upperCI )
        names(paramFixed) = varName

        fitFixed = mle2( minuslogl = negLogL_for_mle2,
                        start = paramGuess,
                        fixed = paramFixed,
                        data  = list( data = data) )

        nllFixed = fitFixed@min
        chiSq    = 2*( nllFixed - nllFull)
        p        = pchisq( chiSq, df = 1)

        upperChiSq   = chiSq
        upperP       = p
        upperConverg = fitFixed@details$convergence 

        if ( verbose ) {
            cat('upperCI = ',upperCI,'\n')
            cat('chiSq   = ',chiSq,  '\n')
            cat('p       = ',p,      '\n')
            cat('delta   = ',delta,  '\n\n')
        }


        # if ( p < alphaLevel && delta > 0 ) { } # keep going

        if ( p > alphaLevel && delta > 0 ) { delta = - delta/2.0 } # switch directions

        # if ( p > alphaLevel && delta < 0 ) { } # keep going 

        if ( p < alphaLevel && delta < 0 ) { delta = - delta/2.0 } # switch directions

    }

    ret = data.frame( name = varName, lowerCI = lowerCI, lowerChiSq = lowerChiSq, lowerP = lowerP, lowerConverg = lowerConverg,
                      upperCI = upperCI, upperChiSq = upperChiSq, upperP = upperP, upperConverg = upperConverg )

    return(ret)
}

########################################################
# Define some auxiliary functions
#
########################################################

record_param_set = function( result, filename)
{
    if ( !file.exists(filename) ) {
        output             = data.frame( matrix(result@fullcoef, nrow=1) ) 
        names(output)      = names(result@fullcoef)
        output$convergence = result@details$convergence
        output$negLogL     = result@min
        output$delta       = 0
        output$chisq       = 0
        output$p           = 0
        output$type        = 'Fit'
        write.csv(output, filename, row.names=F, quote=F)
    } else { 
        output             = data.frame( matrix(result@fullcoef, nrow=1) ) 
        names(output)      = names(result@fullcoef)
        output$convergence = result@details$convergence
        output$negLogL     = result@min
        output$delta       = 0
        output$chisq       = 0
        output$p           = 0
        output$type        = 'Fit'
        write.table(output, filename, row.names=F, append=T, quote=F, col.names=F, sep=',')
    }

}
########################################################
# Perform the optimization
#
########################################################

if ( T ) {
    for ( replicate in 1:10 ) {

        initGuessMin = -1.0
        initGuessMax = 1.0
        paramGuess    = c( beta0 = runif(1, min=initGuessMin, max=initGuessMax ),

                            beta1     = runif(1, min=initGuessMin, max=initGuessMax ),

                            #beta2     = runif(1, min=initGuessMin, max=initGuessMax ),

                            beta4_Col = runif(1, min=initGuessMin, max=initGuessMax ),
                            beta4_Mex = runif(1, min=initGuessMin, max=initGuessMax ),
                            beta4_Swe = runif(1, min=initGuessMin, max=initGuessMax ),
                            beta4_USA = runif(1, min=initGuessMin, max=initGuessMax ),
                            beta4_Ven = runif(1, min=initGuessMin, max=initGuessMax ),

                            #beta4a_NotLA = runif(1, min=initGuessMin, max=initGuessMax ),

                            beta5_Aca = runif(1, min=initGuessMin, max=initGuessMax ),
                            beta5_All = runif(1, min=initGuessMin, max=initGuessMax ),
                            beta5_Cyp = runif(1, min=initGuessMin, max=initGuessMax ),
                            beta5_Gym = runif(1, min=initGuessMin, max=initGuessMax ),
                            beta5_Sil = runif(1, min=initGuessMin, max=initGuessMax ),

                            #beta5a_Specialist = runif(1, min=initGuessMin, max=initGuessMax ),

                            beta6_Mor = runif(1, min=initGuessMin, max=initGuessMax ),
                            beta6_Nom = runif(1, min=initGuessMin, max=initGuessMax ),
                            beta6_Phy = runif(1, min=initGuessMin, max=initGuessMax ),
                            beta6_Rep = runif(1, min=initGuessMin, max=initGuessMax ),

                            beta7_Neo = runif(1, min=initGuessMin, max=initGuessMax ),
                            beta7_Pro = runif(1, min=initGuessMin, max=initGuessMax ),
                            beta7_Tet = runif(1, min=initGuessMin, max=initGuessMax ),

                            #beta8_MOR = runif(1, min=initGuessMin, max=initGuessMax ),
                            #beta8_MOL = runif(1, min=initGuessMin, max=initGuessMax ),

                            beta11_MOR = runif(1, min=initGuessMin, max=initGuessMax ),
                            beta11_MOL = runif(1, min=initGuessMin, max=initGuessMax ),

                            beta14 = 10,

                            #beta14a = runif(1, min=initGuessMin, max=initGuessMax ),
                            #beta14b = runif(1, min=initGuessMin, max=initGuessMax ),
                            #beta14c = runif(1, min=initGuessMin, max=initGuessMax ),
                            #beta14d = runif(1, min=initGuessMin, max=initGuessMax ),
                            #beta14e = runif(1, min=initGuessMin, max=initGuessMax ),
                            #beta14f = runif(1, min=initGuessMin, max=initGuessMax ),
                            #beta14g = runif(1, min=initGuessMin, max=initGuessMax ),

                            beta15 = runif(1, min=initGuessMin, max=initGuessMax ),

                            shape  = shape_map_inv(1.0) )

        cont         = TRUE
        negLogLMin   = negLogL_for_mle2( paramGuess, data )
        iter         = 0
        paramBest    = paramGuess
        paramCurrent = paramGuess

        while ( cont && iter < 100 ) {
            cat('james2.R : Msg : Running optimzaiton\n')
            result      = mle2( minuslogl = negLogL_for_mle2, 
                                start     = paramCurrent, 
                                method    = methodLoc,
                                optimizer = optimizerLoc,
                                data      = list( data = data),
                                control   = list( maxit = maxitLoc ) )

            cat('iter = ',iter,' negLogL = ',result@min,'\n')
            #record_param_set( result, 'james_record.csv')

            if ( negLogLMin - result@min > 0.01 ) {
                negLogLMin = result@min
                paramBest  = result@fullcoef
                resultBest = result
                iter       = 0
            } else { 
                iter = iter + 1
            }

            paramCurrent = paramBest * runif( length(paramBest), 0.8, 1.2)

        }

        cat('james2.R : Msg : Recording data\n')
        record_param_set( resultBest, resultFilename)
    }
}

########################################################
# Compute the profiles and confidence intervals.
#
########################################################
cat('james2.R : Msg : Computing CI\n')

resultCI2 = confint2( resultBest, alphaLevel=0.9667 )

write.csv(resultCI2, ciFilename, row.names=F, quote=F)

########################################################
# Visualize the results
########################################################
#yHat = model(result$par, data)

#plot(Y~X, data=data)
#lines(data$X, yHat, col='red', lwd=3)

