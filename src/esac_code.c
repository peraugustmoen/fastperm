#include "header.h"

void internal_threshold_matrix(double * matrix, int r1, int c1, double a, double nu_a, int previously_tresholded,
                                double prev_nu_a ){

    // should be possible to optimize this code a bit
    double a_sq = a*a;
    double true_val = 0;
    if (previously_tresholded)
    {
        for (int i = 0; i < r1; ++i)
        {
            for (int j = 0; j < c1; ++j)
            {
                if (fabs(matrix[cord_spec(i, j, r1)])>1e-10)
                {
                    true_val = matrix[cord_spec(i, j, r1)]+prev_nu_a;

                    if (true_val>a_sq)
                    {
                        matrix[cord_spec(i, j, r1)] = true_val - nu_a;
                    }
                    else{
                        matrix[cord_spec(i, j, r1)] = 0.0;
                    }
                }

            }

        }
    }
    else{
        for (int i = 0; i < r1; ++i)
        {
            for (int j = 0; j < c1; ++j)
            {
                if (fabs(matrix[cord_spec(i, j, r1)])>a)
                {
                    matrix[cord_spec(i, j, r1)] = matrix[cord_spec(i, j, r1)]*matrix[cord_spec(i, j, r1)] - nu_a;

                }
                else{
                    matrix[cord_spec(i, j, r1)] = 0.0;
                }

            }

        }
    }


}

void internal_colSum(double * matrix, int r1, int c1, double * vector){
    //memset(vector, 0, c1);
    for (int j = 0; j < c1; ++j)
    {
        vector[j]=0.0;
        for (int i = 0; i < r1; ++i)
        {
            vector[j]+=matrix[cord_spec(i, j, r1)];
        }
    }

}
void internal_check_segment(double * cumsums, double * cusum, int * maxpos, double * maximum, int * maxa_pos,
                            int s, int e, int p, double * vector, double * thresholds, double * thresholds_test,
                            double * as, double * nu_as, int len_as, double * tmpvec, int twologn, int * ts,
                            int fast, int debug){
    // require as[0] = 0
    if(debug){
      Rprintf("checking segment (%d,%d]. FAST = %d.\n", s,e, fast);

    }

    *maxpos  = s+(e-s)/2;
    *maximum = NEGINF+1;
    *maxa_pos = 0;
    if(e-s<2){
        return;
    }

    int k = 0;



    //double prev_nu_a = as[0];
    //double a;
    //double nu_a;
    double tmp = 0;
    double tmpmax = 0;
    int tmppos = 0;
    int tmpapos = 0;
    //double a_tmp=-100000;
    //int pos_a_tmp = 0;
    //double * val=0;
    //int n = e-s;

    //int detected = 0;
    double tmp2 =0;
    int localdetected = 0;
    int end = 1;
    int detected = 0;
    int j = 0;
    if(fast){
        k = (s+e)/2;

        // double check this singleCUSUM function
        singleCUSUM(cumsums, cusum, s, e, p, k);
        memset(tmpvec, 0, sizeof(double)*len_as);
        
        // first aggregate thresholded CUSUMs
        j = 0;
        for (int i = 0; i < p; ++i){
            tmp2 = cusum[cord_spec(i,j, p)]*cusum[cord_spec(i,j, p)];
            for (int z = 0; z < len_as; ++z){
                if (fabs(cusum[cord_spec(i,j, p)])>as[z])
                {
                    tmpvec[z] += tmp2 - nu_as[z];
                }
                else{
                    break;
                }
                
            }
        }

        // check if T >0
        for (int z = 0; z < len_as; ++z)
        {
            tmp = tmpvec[z] - thresholds_test[z];
            if(ts[z] <= twologn){
                break;
            }
            else if(tmp>0){
                detected = 1;
                break;
            }
        }

        // if T<= 0, check if P >0:
        if(twologn>0 && detected ==0){
            sort_k_largest_abs(cusum + cord_spec(0,j,p) , twologn, 0, p);
            double cumsum = 0;
            int prev = 0;

            for (int z = (len_as-1); z >= 0; --z)
            {
                if(ts[z] > twologn){
                    break;
                }
     
                if(z==(len_as-1)){
                    prev = -1;
                }
                for (int hh = prev+1; hh < ts[z]; ++hh)
                {
                    cumsum += cusum[cord_spec(hh,j,p)]*cusum[cord_spec(hh,j,p)];
                }
                prev = ts[z]-1;
                tmp = cumsum;

                if( (tmp ) >thresholds_test[z]){
                    detected = 1;
                    break;
                }
            }
        }

        if(detected && debug){
            Rprintf("Checked pos k = %d in [%d, %d) and found chgpt\n", k, s, e);
        }

        if(detected){
            // check all positions
            CUSUM(cumsums, cusum, s, e, p);
            end = e-s-1;
            for (j = 0; j < end; ++j){
                memset(tmpvec, 0, sizeof(double)*len_as);
                
                // first aggregate thresholded CUSUMs
                for (int i = 0; i < p; ++i){
                    tmp2 = cusum[cord_spec(i,j, p)]*cusum[cord_spec(i,j, p)];
                    for (int z = 0; z < len_as; ++z){
                        if (fabs(cusum[cord_spec(i,j, p)])>as[z])
                        {
                            tmpvec[z] += tmp2 - nu_as[z];
                        }
                        else{
                            break;
                        }
                        
                    }
                }

                if(debug){
                    Rprintf("Computing S statistic in k = %d in [%d, %d)\n", s+j+1, s, e);
                }
                tmpmax = NEGINF;
                tmpapos = 0;
                for (int z = 0; z < len_as; ++z)
                {
                    /*if(fabs(tmpvec[z])<1e-10){
                        continue;
                    }*/
                    
                    tmp = tmpvec[z] - thresholds[z];
                    if(tmp > tmpmax){
                        tmpmax = tmp;
                        tmpapos = z;

                    }
                    if (tmp> *maximum)
                    {
                        *maximum= tmp;
                        *maxpos = s+j+1;
                        *maxa_pos = z;
                    }
                    
                }
                if(debug){
                    Rprintf("S stat = %f at apos %d\n", tmpmax, tmpapos);
                }
            }
            if(debug){
                Rprintf("for segment (%d,%d] we found maximum in %d with val %f\n", s,e,*maxpos, *maximum);
            }
        }

        return;
        
    }



    // this runs only if fast = FALSE
    

    if(debug){
        Rprintf("Checking all poses in [%d, %d)\n", s, e);
    }

    // now full thing. 
    CUSUM(cumsums, cusum, s, e, p);
    end = e-s-1;
    for (j = 0; j < end; ++j){
        memset(tmpvec, 0, sizeof(double)*len_as);
        
        // first aggregate thresholded CUSUMs
        for (int i = 0; i < p; ++i){
            tmp2 = cusum[cord_spec(i,j, p)]*cusum[cord_spec(i,j, p)];
            for (int z = 0; z < len_as; ++z){
                if (fabs(cusum[cord_spec(i,j, p)])>as[z])
                {
                    tmpvec[z] += tmp2 - nu_as[z];
                }
                else{
                    break;
                }
                
            }
        }

        if(debug){
            Rprintf("Computing S statistic in k = %d in [%d, %d)\n", s+j+1, s, e);
        }
        tmpmax = NEGINF;
        tmpapos = 0;
        for (int z = 0; z < len_as; ++z)
        {
            if(fabs(tmpvec[z])<1e-10){
                continue;
            }
            //tmp = tmpvec[z] - thresholds[z];
            //if(tmpvec[z] > thresholds_test[z]){
                tmp = tmpvec[z] - thresholds[z];
                if(tmp > tmpmax){
                    tmpmax = tmp;
                    tmpapos = z;

                }
                if (tmp> *maximum)
                {
                    *maximum= tmp;
                    *maxpos = s+j+1;
                    *maxa_pos = z;
                }
            //}
            
        }
        if(debug){
            Rprintf("S stat = %f at apos %d\n", tmpmax, tmpapos);
        }

        // check if T >0
        localdetected = 0;

        if(!detected){
            for (int z = 0; z < len_as; ++z)
            {
                /*tmp = tmpvec[z] - thresholds[z];
                if (tmp> *maximum)
                {
                    *maximum= tmp;
                    *maxpos = s+j+1;
                    *maxa_pos = z;
                }*/
                tmp = tmpvec[z] - thresholds_test[z];
                if(tmp>0 && ts[z] > twologn){
                    localdetected = 1;
                }
            }

            // if T<= 0, check if P >0:
            if(twologn>0 && localdetected ==0){
                sort_k_largest_abs(cusum + cord_spec(0,j,p) , twologn, 0, p);
                //partial_quicksort(cusum + cord_spec(0,j,p) , p , twologn);
                double cumsum = 0;
                int prev = 0;

                for (int z = (len_as-1); z >= 0; --z)
                {
                    if(ts[z] > twologn){
                        break;
                    }
         
                    if(z==(len_as-1)){
                        prev = -1;
                    }
                    for (int hh = prev+1; hh < ts[z]; ++hh)
                    {
                        cumsum += cusum[cord_spec(hh,j,p)]*cusum[cord_spec(hh,j,p)];
                    }
                    prev = ts[z]-1;
                    tmp = cumsum;

                    if(tmp >thresholds_test[z]){
                        localdetected = 1;
                    }
                }
            }
        
            if(localdetected){
                detected=1;
                if (debug)
                {
                    Rprintf("detected chgpt at pos %d in (%d,%d]\n", s+j+1,s,e);
                }
            }
        }
    }

    if(detected==0){
        *maximum = NEGINF+1;
    }
    if(debug){
      Rprintf("for segment (%d,%d] we found maximum in %d with val %f\n", s,e,*maxpos, *maximum);
    }

    return;
}

void internal_find_changepoint(double * cumsums, double * cusum, int * maxpos, double * maximum, int * maxa_pos,
                            int s, int e, int p, double * vector, double * thresholds, double * as, double * nu_as, int len_as, double * tmpvec, int twologn, int * ts,
                            int debug){
    // require as[0] = 0
    /*if(debug){
      Rprintf("checking segment (%d,%d]\n", s,e);

    }*/

    *maxpos  = s+(e-s)/2;
    *maximum = NEGINF+1;
    *maxa_pos = 0;
    if(e-s<2){
        return;
    }


    //double prev_nu_a = as[0];
    //double a;
    //double nu_a;
    double tmp = 0;
    //double a_tmp=-100000;
    //int pos_a_tmp = 0;
    //double * val=0;
    //int n = e-s;

    //int detected = 0;
    double tmp2 =0;
    int end = 1;

    CUSUM(cumsums, cusum, s, e, p);
    end = e-s-1;
    for (int j = 0; j < end; ++j){
        memset(tmpvec, 0, sizeof(double)*len_as);
        
        // first aggregate thresholded CUSUMs
        for (int i = 0; i < p; ++i){
            tmp2 = cusum[cord_spec(i,j, p)]*cusum[cord_spec(i,j, p)];
            for (int z = 0; z < len_as; ++z){
                if (fabs(cusum[cord_spec(i,j, p)])>as[z])
                {
                    tmpvec[z] += tmp2 - nu_as[z];
                }
                else{
                    break;
                }
                
            }
        }


        if(debug){
            Rprintf("Checked pos k = %d in [%d, %d) and found chgpt", s+j+1, s, e);
        }
        for (int z = 0; z < len_as; ++z)
        {
            if(fabs(tmpvec[z])<1e-10){
                continue;
            }
            //tmp = tmpvec[z] - thresholds[z];
            //if(tmpvec[z] > thresholds_test[z]){
                tmp = tmpvec[z] - thresholds[z];
                if (tmp> *maximum)
                {
                    *maximum= tmp;
                    *maxpos = s+j+1;
                    *maxa_pos = z;
                }
            //}
            
        }
    }

    return;
}

/*
void internal_find_maxpos(double * cumsums, double * cusum, int * maxpos, double * maximum, int * maxa_pos,
                            int s, int e, int jval, int * k_max, int * segstarts, 
                            int * lens, int lenLens, int p, int n, double * vector, double * thresholds,
                            double * as, double * nu_as, int len_as, double * tmpvec,  int * ts,
                            int * maxposes, double * maxvalues,  int * maxas, int debug){
    // require as[0] = 0
    if(debug){
      Rprintf("checking segment (%d,%d]\n", s,e);

    }

    *maxpos  = s+(e-s)/2;
    *maximum = NEGINF+1;
    *maxa_pos = 0;
    *k_max = 0;
    double tmpmax;
    int tmppos;
    int tmpmaxapos;

    int ss; 
    int ee;
    int i;
    int j = jval;
    int len = lens[j];
    int k;
    int tmp2;
    int tmp;

    if(e-s<2){
        return;
    }

    for (k = 0; k < n; ++k)
    {

        i = segstarts[cord_spec(k,j,n)];
        if(debug){
            Rprintf("i= %d\n", i);
        }
        if(i>e-len || i<-1){
            if(debug){
                Rprintf("i= %d is skipped\n", i);
            }
            break;
        }
        else if (i<s)
        {
            continue;
        }

        if(debug){
            //Rprintf("maxvalues[%d, %d] = %f\n", k , j , maxvalues[cord_spec(k,j,n)]);
        }
        


        if((NEGINF +1) <= maxvalues[cord_spec(k,j,n)] && maxvalues[cord_spec(k,j,n)] <=(NEGINF+2)){
            ss = i;
            ee = i + len;
            tmpmax = NEGINF+3;
            tmppos =  (ss+ee)/2;
            tmpmaxapos = 1;
            CUSUM(cumsums, cusum, ss, ee, p);
             if(debug){
                Rprintf("segment (%d,%d] (k=%d, j=%d) not inspected, now checking!\n", i, i+len,k,j);
             }
            for (int j = 0; j <ee-ss-1; ++j){
                memset(tmpvec, 0, sizeof(double)*len_as);
                
                // first aggregate thresholded CUSUMs
                for (int i = 0; i < p; ++i){
                    tmp2 = cusum[cord_spec(i,j, p)]*cusum[cord_spec(i,j, p)];
                    for (int z = 0; z < len_as; ++z){
                        if (fabs(cusum[cord_spec(i,j, p)])>as[z])
                        {
                            tmpvec[z] += tmp2 - nu_as[z];
                        }
                        else{
                            break;
                        }
                        
                    }
                }
                for (int z = 0; z < len_as; ++z)
                {
                    if(fabs(tmpvec[z])<1e-10){
                        continue;
                    }
                    //tmp = tmpvec[z] - thresholds[z];
                    //if(tmpvec[z] > thresholds_test[z]){
                    tmp = tmpvec[z] - thresholds[z];
                    if (tmp> tmpmax)
                    {
                        tmpmax= tmp;
                        tmppos = s+j+1;
                        tmpmaxapos = z;
                    }

                        
                }
            }

            maxvalues[cord_spec(k,j,n)] = tmpmax;
            maxpos[cord_spec(k,j,n)] = tmppos;
            maxas[cord_spec(k,j,n)] = tmpmaxapos;
        }

        else if(debug && maxvalues[cord_spec(k,j,n)] > (NEGINF+2)){
            Rprintf("segment (%d,%d] (k=%d, j=%d) already inspected, with max val %f in %d\n", i, i+len,k,j,
            maxvalues[cord_spec(k,j,n)],maxpos[cord_spec(k,j,n)]);
        }
        else if(maxvalues[cord_spec(k,j,n)] <=(NEGINF+1)){
            continue;
        }
        tmp = maxvalues[cord_spec(k,j,n)];

        if(tmp>*maximum){
            *maximum = tmp;
            *maxpos = maxpos[cord_spec(k,j,n)];
            *maxa_pos = maxas[cord_spec(k,j,n)];
            *k_max=k;
            //found=1;
        }


    }

    
    if(debug){
      Rprintf("for segment (%d,%d] we found maximum in %d with val %f\n", s,e,*maxpos, *maximum);
    }

    return;
}*/


void cESAC_call(double * x, int s, int e, int n, int p, int depth, int* changepoints,int* changepoint_counter_ptr, int * depthcounter,
                double * thresholds, double * thresholds_test, double *cumsums, int* lens, int lenLens, double * as, double * nu_as, int len_as,
                int * segstarts, double * maxvalues, int* maxpos,int * maxas, int K, double * cusum,
                double * vector, int * coordchg, double * maxval, int * startpoints, int* endpoints, int* maxaposes,
                double * tmpvec, int twologn, int * ts, int fast, int trim, int NOT,int midpoint, int debug){
    if(debug){
        Rprintf("cESAC_call! s=%d, e=%d\n", s, e);
    }

    if(e-s<2*lens[0]){
        //Rprintf("segment too short\n");
        return;
    }
    int argmax = s+1;
    double maximum = NEGINF;
    int maxa_pos = 0;

    //int tmpargmax;
    //double tmpmaximum;
    double tmp;
    int len;
    //int found=0;

    int i;

    int j_max=0;
    int k_max = 0;
    int k;
    int j;

    for ( j = 0; j < lenLens; ++j)
    {
        len = lens[j];
        if(debug){
            Rprintf("j=%d, len = %d\n", j, len);
        }



        if(e-s<2*len){
            break;
        }

        for (k = 0; k < n; ++k)
        {
            i = segstarts[cord_spec(k,j,n)];
            if(debug){
                Rprintf("i= %d\n", i);
            }
            if(i>e-2*len || i<-1){
                if(debug){
                    Rprintf("i= %d is skipped\n", i);
                }
                break;
            }
            else if (i<s)
            {
                continue;
            }

            if(debug){
                //Rprintf("maxvalues[%d, %d] = %f\n", k , j , maxvalues[cord_spec(k,j,n)]);
            }

            if(maxvalues[cord_spec(k,j,n)]<=NEGINF){
                 if(debug){
                    Rprintf("segment (%d,%d] (k=%d, j=%d) not inspected, now checking!\n", i, i+2*len,k,j);
                 }
                //this segment not computed!
                //inspectOnSegment(cumsums, cusum, &tmpargmax, &tmpmaximum, s, e, p, lambda,
                 //   eps, maxiter, projvec, cusum_proj);
                 //double * cumsums, double * cusum, int * maxpos, double * maximum,
                 internal_check_segment(cumsums, cusum, &(maxpos[cord_spec(k,j,n)]), &(maxvalues[cord_spec(k,j,n)]), &(maxas[cord_spec(k,j,n)]),
                            i, i+2*len, p, vector, thresholds, thresholds_test, as, nu_as,len_as, tmpvec, twologn, ts,fast,debug);
                 //internal_inspectOnSegment(cumsums, cusum, &(maxpos[cord_spec(k,j,n)]), &(maxvalues[cord_spec(k,j,n)]), i, i+len, p,
                 //lambda,
                //      eps, maxiter, mhat, mhatprod, v, v2,debug);
            }
            else if(debug){
                Rprintf("segment (%d,%d] (k=%d, j=%d) already inspected, with max val %f in %d\n", i, i+2*len,k,j,
                maxvalues[cord_spec(k,j,n)],maxpos[cord_spec(k,j,n)]);
            }
            tmp = maxvalues[cord_spec(k,j,n)];

            if(tmp>maximum){
                maximum = tmp;
                argmax = maxpos[cord_spec(k,j,n)];
                maxa_pos = maxas[cord_spec(k,j,n)];
                j_max = j;
                k_max=k;
                //found=1;
            }


        }

        if( NOT && maximum>NEGINF+1){
            break;
        }
/*        if(maximum>0.0){
            break;
        }
*/    }
    if(debug){
        Rprintf("maximum=%f\n", maximum);
    }

    /*if(fast && maximum > NEGINF+1){
        // find maximizing pos
    }*/
    if(maximum >NEGINF+1){
        if(debug){
          Rprintf("!!!!!! declared change-point in %d. val = %f. (s,e] = (%d,%d]\n", argmax, maximum,
          segstarts[cord_spec(k_max,j_max,n)], 2*lens[j_max]+ segstarts[cord_spec(k_max,j_max,n)]);
          Rprintf("changeptcounter = %d\n", *changepoint_counter_ptr);
        }


        // identify in which coordinates the change happens:
        //i = segstarts[cord_spec(k_max,j_max,n)];
        //len = lens[j_max];
        //int ss = i;
        //int ee = i+len;
        if(maxa_pos==0){
            for (int zz = 0; zz < p; ++zz)
            {
                coordchg[cord_spec(zz,*changepoint_counter_ptr, p)]=1;
            }
        }
        else{
            i = segstarts[cord_spec(k_max,j_max,n)];
            len = lens[j_max];
            int ss = i;
            int ee = i+2*len;
            //CUSUM(cumsums, cusum, ss, ee, p);
            singleCUSUM(cumsums, cusum, ss, ee, p, argmax);
            //internal_threshold_matrix(&(cusum[cord_spec(0,argmax-ss-1,p)]), p, 1, as[maxa_pos],  nu_as[maxa_pos], 0,
            //                        0 );
            internal_threshold_matrix(cusum, p, 1, as[maxa_pos],  nu_as[maxa_pos], 0,
                                    0 );


            for (int zz = 0; zz < p; ++zz)
            {
                if(cusum[cord_spec(0,argmax-ss-1,p)+zz]>1e-10){
                    coordchg[cord_spec(zz,*changepoint_counter_ptr, p)]=1;
                }
            }

        }

        changepoints[*changepoint_counter_ptr] = argmax;
        depthcounter[*changepoint_counter_ptr] = depth;
        maxval[*changepoint_counter_ptr] = maximum;
        startpoints[*changepoint_counter_ptr] = segstarts[cord_spec(k_max,j_max,n)];
        endpoints[*changepoint_counter_ptr] = segstarts[cord_spec(k_max,j_max,n)] + 2*lens[j_max];
        maxaposes[*changepoint_counter_ptr] = maxa_pos;
        if(midpoint){
            changepoints[*changepoint_counter_ptr] = (startpoints[*changepoint_counter_ptr] + endpoints[*changepoint_counter_ptr])/2;
            argmax = changepoints[*changepoint_counter_ptr];
        }
        
        //int cut = (int) ( lens[j_max] * 2*cutoff) ;
        //Rprintf("split at %d. s = %d, e = %d. cut  = %d\n", argmax, s, e,cut);
        //cESAC_call(double * x, int s, int e, int n, int p, int depth, int* changepoints,int* changepoint_counter_ptr, int * depthcounter,
        //        double threshold_d, double threshold_s , double *cumsums, int* lens, int lenLens, double * as, double * nu_as, int len_as,
        //        int * segstarts, double * maxvalues, int* maxpos, int K, double * cusum, double * matrix,
        //        double * vector, int * coordchg, int debug)

        int s1 = s;
        int e1 = argmax;
        int s2 = argmax;
        int e2 = e;

        if(trim){
            e1 = startpoints[*changepoint_counter_ptr]+1;
            s2 = endpoints[*changepoint_counter_ptr]-1;
        }

        (*changepoint_counter_ptr)++;
        if(*changepoint_counter_ptr >n){
          return;
        }

        cESAC_call(x, s1, e1, n, p, depth+1, changepoints,changepoint_counter_ptr,  depthcounter,
                thresholds, thresholds_test, cumsums, lens, lenLens,  as, nu_as, len_as,
                segstarts, maxvalues, maxpos,maxas, K, cusum, vector, coordchg, maxval, 
                startpoints, endpoints, maxaposes, tmpvec,twologn, ts,fast, trim, NOT,midpoint,debug);
        cESAC_call(x, s2, e2,n, p, depth+1, changepoints,changepoint_counter_ptr,  depthcounter,
                thresholds , thresholds_test, cumsums, lens, lenLens,  as, nu_as, len_as,
                segstarts, maxvalues, maxpos,maxas, K, cusum, vector, coordchg, maxval,
                startpoints, endpoints, maxaposes, tmpvec, twologn, ts,fast, trim, NOT,midpoint, debug);

    }

    return;
}

SEXP cESAC(SEXP XI,SEXP nI, SEXP pI,SEXP thresholdsI, SEXP thresholds_testI, SEXP lensI,SEXP lenLensI,SEXP KI,
    SEXP asI, SEXP nu_asI, SEXP len_asI, SEXP twolognI, SEXP tsI, SEXP fastI, SEXP rescale_variance_boolI, 
    SEXP trimI, SEXP NOTI,SEXP midpointI,SEXP debugI){
    // X : p \times n
    PROTECT(XI);
    PROTECT(thresholdsI);
    PROTECT(thresholds_testI);
    PROTECT(lensI);
    PROTECT(asI);
    PROTECT(nu_asI);
    PROTECT(nI);
    PROTECT(pI);
    PROTECT(lenLensI);
    PROTECT(KI);
    PROTECT(len_asI);
    PROTECT(twolognI);
    PROTECT(tsI);
    PROTECT(debugI);
    PROTECT(rescale_variance_boolI);
    PROTECT(trimI);
    PROTECT(fastI);
    PROTECT(NOTI);
    PROTECT(midpointI);

    double * X = REAL(XI);
    int n = *(INTEGER(nI));
    int p = *(INTEGER(pI));
    double *thresholds = (REAL(thresholdsI));
    double *thresholds_test = (REAL(thresholds_testI));
    int *lens = INTEGER(lensI); /////// dobbeltsjekk at int er like stor som INTEGER!!!!!!!
    int lenLens = *(INTEGER(lenLensI));
    int K = *(INTEGER(KI));
    double * as = REAL(asI);
    double * nu_as = REAL(nu_asI);
    int len_as = *(INTEGER(len_asI));
    int twologn = * (INTEGER(twolognI));
    int * ts = INTEGER(tsI);
    int debug = *INTEGER(debugI);
    int fast = *INTEGER(fastI);
    int rescale_variance_bool = *INTEGER(rescale_variance_boolI);
    int trim = *INTEGER(trimI);
    int NOT = *INTEGER(NOTI);
    int midpoint = *INTEGER(midpointI);

    if(debug){
      Rprintf("p = %d\n", p);
      Rprintf("n = %d\n", n);
    }
    SEXP out = PROTECT(allocVector(INTSXP, n));
    int * changepoints = INTEGER(out); //pointer to array
    //memset(changepoints, -1, sizeof(int)*n);
    for (size_t i = 0; i < n; i++) {
      changepoints[i] = -1;
    }
    int changepoint_counter = 0;
    int* changepoint_counter_ptr = &changepoint_counter;
    SEXP maxvalSEXP = PROTECT(allocVector(REALSXP, n));
    double * maxval = REAL(maxvalSEXP); //pointer to array
    memset(maxval, 0, sizeof(double)*n);
    SEXP startpointsSEXP = PROTECT(allocVector(INTSXP, n));
    int * startpoints = INTEGER(startpointsSEXP); //pointer to array
    memset(startpoints, 0, sizeof(int)*n);
    SEXP endpointsSEXP = PROTECT(allocVector(INTSXP, n));
    int * endpoints = INTEGER(endpointsSEXP); //pointer to array
    memset(endpoints, 0, sizeof(int)*n);
    SEXP maxaposesSEXP = PROTECT(allocVector(INTSXP, n));
    int * maxaposes = INTEGER(maxaposesSEXP); //pointer to array
    memset(maxaposes, 0, sizeof(int)*n);
    SEXP depthcounterSEXP= PROTECT(allocVector(INTSXP, n));
    int * depthcounter = INTEGER(depthcounterSEXP); //pointer to array
    memset(depthcounter, 0, sizeof(int)*n);
    SEXP coordschgSEXP = PROTECT(allocVector(INTSXP, n*p));
    int * coordchg = INTEGER(coordschgSEXP);
    memset(coordchg, 0, sizeof(int)*n*p);
    // first we compute all the cumulative sums of all
    // coordinates:

    int maxlen = p;
    int minlen = n;
    if(n>p){
        maxlen = n;
        minlen = p;
    }
    SEXP vectorSEXP = PROTECT(allocVector(REALSXP, maxlen));


    double * vector = REAL(vectorSEXP);
    memset(vector, 0, sizeof(double)*maxlen);

    SEXP scalesSEXP = PROTECT(allocVector(REALSXP, p));
    double * scales = REAL(scalesSEXP);
    for (int i = 0; i < p; ++i)
    {
        scales[i] = 1;
    }

    if (rescale_variance_bool)
    {
        rescale_variance(X, scales, n, p, vector);
    }

    SEXP cumsumsSEXP = PROTECT(allocVector(REALSXP, p * (n+1)));
    double *cumsums = REAL(cumsumsSEXP); // p \times (n+1). first col is 0
    memset(cumsums, 0, p*(n+1)*sizeof(double));

    for (int j = 1; j <=n; ++j)
    {
        for (int i = 0; i < p; ++i)
        {
            //#define cord_spec(r,c, D) ((r) + (D)*(c))
            cumsums[cord_spec(i,j,p)] = X[cord_spec(i,j-1, p)] +cumsums[cord_spec(i,j-1, p)];
        }
    }

    SEXP cusumSEXP = PROTECT(allocVector(REALSXP, p*(n)));
    double * cusum = REAL(cusumSEXP);
    memset(cusum, 0, sizeof(double)*p*(n));
    




    // record all segment starts
    // we don't compute all cusums here but do it on the fly so we dont have to do anything
    // unecessary


    SEXP maxvaluesSEXP = PROTECT(allocVector(REALSXP, n*lenLens));
    double * maxvalues = REAL(maxvaluesSEXP);
    //memset(maxvalues, 0, sizeof(double)*n*lenLens);
    for (int i = 0; i < n*lenLens; ++i)
    {
        maxvalues[i] = NEGINF;
    }
    SEXP tmpvecSEXP = PROTECT(allocVector(REALSXP, len_as));
    double * tmpvec = REAL(tmpvecSEXP);
    memset(tmpvec, 0, sizeof(double)*len_as);
    SEXP maxposSEXP = PROTECT(allocVector(INTSXP, n*lenLens));
    int * maxpos = INTEGER(maxposSEXP);
    memset(maxpos, 0, sizeof(int)*n*lenLens);
    SEXP maxasSEXP = PROTECT(allocVector(INTSXP, n*lenLens));
    int * maxas = INTEGER(maxasSEXP);
    memset(maxas, 0, sizeof(int)*n*lenLens);
    SEXP segstartsSEXP = PROTECT(allocVector(INTSXP, n*lenLens));
    int * segstarts = INTEGER(segstartsSEXP);
    //memset(segstarts, -2, sizeof(int)*n*lenLens);
    for (size_t i = 0; i < n*lenLens; i++) {
      segstarts[i] = -2;
    }
    int len;
    int jump;
    int counter = 0;
    for (int j = 0; j < lenLens; ++j)
    {
        counter=0;

        len = lens[j]; //len here is -1 of the len in text
        jump = len/K;
        if(jump<1){
            jump=1;
        }

        for (int i = -1; i < (n-2*len); i+=jump)
        {
            //cord_spec(r,c, D) ((r) + (D)*(c))
            //compute_cusum(cumsum, i, i+len, &(maxpos[cord_spec(i,j,n)]), &(maxvalues[cord_spec(i,j,n)]));
            segstarts[cord_spec(counter++,j,n)] = i;
            if(debug){
                //Rprintf("segstarts[%d, %d] = %d\n",counter-1, j, i );
            }


        }
    }

    if(debug){
      for (int j = 0; j < lenLens; ++j)
      {
          //counter=0;



          for (int i = 0; i < n; i++)
          {
              //cord_spec(r,c, D) ((r) + (D)*(c))
              //compute_cusum(cumsum, i, i+len, &(maxpos[cord_spec(i,j,n)]), &(maxvalues[cord_spec(i,j,n)]));
              //segstarts[cord_spec(counter++,j,n)] = i;
              if(debug){
                  Rprintf("segstarts[%d, %d] = %d\n",i, j, segstarts[cord_spec(i,j,n)] );
              }


          }
      }
    }




    cESAC_call(X, -1, n-1, n, p, 1, changepoints,changepoint_counter_ptr,  depthcounter,
                thresholds , thresholds_test, cumsums, lens, lenLens,  as, nu_as, len_as,
                segstarts, maxvalues, maxpos, maxas, K, cusum, vector, coordchg,maxval, 
                startpoints, endpoints, maxaposes, tmpvec, twologn, ts,fast,trim,NOT, midpoint,debug);

/*    cInspect_call(X, -1, n-1, n, p, 1, changepoints, changepoint_counter_ptr, depthcounter,
                maxval, xi, cumsums, lens, lenLens, lambda,
                eps, maxiter, segstarts, maxvalues, maxpos, K, cusum, mhat,
                mhatprod, v, v2,debug,coordchg);
*/
    // return:
    //Rprintf("found %d changepoints", *changepoint_counter_ptr);
    SEXP changepointnumSEXP = PROTECT(allocVector(INTSXP, 1));
    int * changepointnum = INTEGER(changepointnumSEXP);
    *changepointnum = changepoint_counter;

    SEXP ret = PROTECT(allocVector(VECSXP, 9)); //the list to be returned in R
    SET_VECTOR_ELT(ret, 0, out);
    SET_VECTOR_ELT(ret, 1, maxvalSEXP);
    SET_VECTOR_ELT(ret, 2, depthcounterSEXP);
    SET_VECTOR_ELT(ret, 3, coordschgSEXP);
    SET_VECTOR_ELT(ret, 4, startpointsSEXP);
    SET_VECTOR_ELT(ret, 5, endpointsSEXP);
    SET_VECTOR_ELT(ret, 6, maxaposesSEXP);
    SET_VECTOR_ELT(ret, 7, changepointnumSEXP);
    SET_VECTOR_ELT(ret, 8, scalesSEXP);

    // creating list of names/titles to be returned in the output list
    SEXP names = PROTECT(allocVector(STRSXP, 9));
    //SET_STRING_ELT(names, 0, mkChar("CUSUM"));
    SET_STRING_ELT(names, 0, mkChar("changepoints"));
    SET_STRING_ELT(names, 1, mkChar("CUSUMval"));
    SET_STRING_ELT(names, 2, mkChar("depth"));
    SET_STRING_ELT(names, 3, mkChar("coordinate"));
    SET_STRING_ELT(names, 4, mkChar("startpoints"));
    SET_STRING_ELT(names, 5, mkChar("endpoints"));
    SET_STRING_ELT(names, 6, mkChar("maxaposes"));
    SET_STRING_ELT(names, 7, mkChar("changepointnumber"));
    SET_STRING_ELT(names, 8, mkChar("scales"));

    setAttrib(ret, R_NamesSymbol, names);

    UNPROTECT(38);
    return ret;
}


SEXP cESAC_calibrate(SEXP nI, SEXP pI, SEXP NI, SEXP tolnI,SEXP lensI,SEXP lenLensI,SEXP KI,
    SEXP asI, SEXP nu_asI, SEXP len_asI, SEXP twolognI, SEXP tsI, SEXP fastI, SEXP rescale_variance_boolI,
    SEXP debugI){
    // X : p \times n
    PROTECT(lensI);
    PROTECT(asI);
    PROTECT(nu_asI);
    PROTECT(nI);
    PROTECT(pI);
    PROTECT(NI);
    PROTECT(tolnI);
    PROTECT(lenLensI);
    PROTECT(KI);
    PROTECT(len_asI);
    PROTECT(twolognI);
    PROTECT(tsI);
    PROTECT(rescale_variance_boolI);
    PROTECT(debugI);

    PROTECT(fastI);

    int n = *(INTEGER(nI));
    int p = *(INTEGER(pI));
    int N = *(INTEGER(NI));
    int toln = *(INTEGER(tolnI));
    int *lens = INTEGER(lensI); /////// dobbeltsjekk at int er like stor som INTEGER!!!!!!!
    int lenLens = *(INTEGER(lenLensI));
    int K = *(INTEGER(KI));
    double * as = REAL(asI);
    double * nu_as = REAL(nu_asI);
    int len_as = *(INTEGER(len_asI));
    int twologn = * (INTEGER(twolognI));
    int * ts = INTEGER(tsI);
    int debug = *INTEGER(debugI);
    int fast = *INTEGER(fastI);
    int rescale_variance_bool = *INTEGER(rescale_variance_boolI);

    if(debug){
      Rprintf("p = %d\n", p);
      Rprintf("n = %d\n", n);
      Rprintf("N = %d\n", N);
      Rprintf("tol*N = %d\n", toln);
    }

    // max vals for no partial sum
    SEXP maxvals1SEXP = PROTECT(allocVector(REALSXP, len_as * N));
    double * maxvals1 = REAL(maxvals1SEXP);
    //memset(maxvals1, 0, sizeof(double)*len_as*N);

    // max vals for partial sum
    SEXP maxvals2SEXP = PROTECT(allocVector(REALSXP, len_as * N));
    double * maxvals2 = REAL(maxvals2SEXP);

    for (int i = 0; i < len_as * N; ++i)
    {
        maxvals2[i] = NEGINF;
        maxvals1[i] = NEGINF;
    }



    SEXP cumsumsSEXP = PROTECT(allocVector(REALSXP, p * (n+1)));
    double *cumsums = REAL(cumsumsSEXP); // p \times (n+1). first col is 0
    memset(cumsums, 0, p*(n+1)*sizeof(double));

 
    SEXP cusumSEXP = PROTECT(allocVector(REALSXP, p*(n)));
    double * cusum = REAL(cusumSEXP);
    memset(cusum, 0, sizeof(double)*p*(n));
    int maxlen = p;
    int minlen = n;
    if(n>p){
        maxlen = n;
        minlen = p;
    }
    SEXP vectorSEXP = PROTECT(allocVector(REALSXP, maxlen));

   
    double * vector = REAL(vectorSEXP);
    memset(vector, 0, sizeof(double)*maxlen);
    SEXP tmpvecSEXP = PROTECT(allocVector(REALSXP, len_as));
    double * tmpvec = REAL(tmpvecSEXP);
    memset(tmpvec, 0, sizeof(double)*len_as);
    SEXP segstartsSEXP = PROTECT(allocVector(INTSXP, n*lenLens));
    int * segstarts = INTEGER(segstartsSEXP);
    //memset(segstarts, -2, sizeof(int)*n*lenLens);
    for (size_t i = 0; i < n*lenLens; i++) {
      segstarts[i] = -2;
    }
    int len;
    int jump;
    int counter = 0;
    for (int j = 0; j < lenLens; ++j)
    {
        counter=0;

        len = lens[j]; //len here is -1 of the len in text
        jump = len/K;
        if(jump<1){
            jump=1;
        }

        for (int i = -1; i < (n-2*len); i+=jump)
        {
            //cord_spec(r,c, D) ((r) + (D)*(c))
            //compute_cusum(cumsum, i, i+len, &(maxpos[cord_spec(i,j,n)]), &(maxvalues[cord_spec(i,j,n)]));
            segstarts[cord_spec(counter++,j,n)] = i;
            /*if(debug){
                //Rprintf("segstarts[%d, %d] = %d\n",counter-1, j, i );
            }*/


        }
    }

    int s = -1;
    int e = n-1;
    int ss;
    int ee;
    int v;
    int ll; 
    int hh;
    int i;
    int j;
    int k;
    int z;
    double tmp2;
    double tmp;

    /*if(debug){
        return NULL;
    }*/
    

    SEXP XSEXP = PROTECT(allocVector(REALSXP, p * (n)));
    double *X = REAL(XSEXP); // p \times (n+1). first col is 0
    memset(X, 0, p*(n)*sizeof(double));

    //SEXP X2 = NULL;

    for (k = 0; k < N; ++k)
    {
        
        // generate X
        GetRNGstate();
        for (j = 0; j < n; ++j)
        {
            for (i = 0; i < p; ++i)
            {
                X[cord_spec(i,j,p)] = norm_rand();
                //X[cord_spec(i,j,p)] = qt(unif_rand(), 10,1,0);
            }
        }
        PutRNGstate();

        // normalize / rescale variance
        if (rescale_variance_bool)
        {
            rescale_variance(X,NULL, n, p, vector);
        }
        
        memset(cumsums, 0, p*sizeof(double));

        for (i = 0; i < p; ++i)
        {
            for (j = 0; j < n; ++j)
            {
                cumsums[cord_spec(i,j+1,p)] = X[cord_spec(i,j,p)] +cumsums[cord_spec(i,j, p)];
            }
        }


        for ( ll = 0; ll < lenLens; ++ll)
        {
            len = lens[ll];
            /*if(debug){
                Rprintf("ll=%d, len = %d\n", ll, len);
            }*/



            if(e-s<2*len){
                break;
            }

            for (hh = 0; hh < n; ++hh)
            {
                ss = segstarts[cord_spec(hh,ll,n)];
                ee = ss+2*len;
                /*if(debug){
                    Rprintf("i= %d\n", hh);
                }*/
                if( ee >n || ss<-1){
                    /*if(debug){
                        Rprintf("i= %d is shhipped\n", hh);
                    }*/
                    break;
                }

                if(fast){
                    v = (ss+ee)/2;
                    singleCUSUM(cumsums, cusum, ss, ee, p, v);
                    memset(tmpvec, 0, sizeof(double)*len_as);
                    
                    // first aggregate thresholded CUSUMs
                    for (i = 0; i < p; ++i){
                        tmp2 = cusum[cord_spec(i,0, p)]*cusum[cord_spec(i,0, p)];
                        for (z = 0; z < len_as; ++z){
                            if (fabs(cusum[cord_spec(i,0, p)])>as[z])
                            {
                                tmpvec[z] += tmp2 - nu_as[z];
                            }
                            else{
                                break;
                            }
                            
                        }
                    }

                    for (int z = 0; z < len_as; ++z)
                    {
                        tmp = tmpvec[z];
                        if(tmp>maxvals1[cord_spec(k, z, N)] ){
                            maxvals1[cord_spec(k, z, N)] = tmp;
                            if(ts[z] > twologn){
                                maxvals2[cord_spec(k,z,N)] = tmp;
                            }
                        }
                    }

                    sort_k_largest_abs(cusum + cord_spec(0,0,p) , twologn, 0, p);
                    //partial_quicksort(cusum + cord_spec(0,j,p) , p , twologn);
                    double cumsum = 0;
                    int prev = 0;

                    for (int z = (len_as-1); z >= 0; --z)
                    {
                        if(ts[z] > twologn){
                            break;
                        }
             
                        if(z==(len_as-1)){
                            prev = -1;
                        }
                        for (int hh = prev+1; hh < ts[z]; ++hh)
                        {
                            cumsum += cusum[cord_spec(hh,0,p)]*cusum[cord_spec(hh,0,p)];
                        }
                        prev = ts[z]-1;
                        tmp = cumsum;

                        if( (tmp ) >maxvals2[cord_spec(k,z,N)]){
                            maxvals2[cord_spec(k,z,N)] = tmp;
                        }
                    }
                    
                }
                else{
                    CUSUM(cumsums, cusum, ss, ee, p);
                    for (v = 0; v < ee-ss-1; ++v)
                    {
                        
                        memset(tmpvec, 0, sizeof(double)*len_as);
                        
                        // first aggregate thresholded CUSUMs
                        for (i = 0; i < p; ++i){
                            tmp2 = cusum[cord_spec(i,v, p)]*cusum[cord_spec(i,v, p)];
                            for (z = 0; z < len_as; ++z){
                                if (fabs(cusum[cord_spec(i,v, p)])>as[z])
                                {
                                    tmpvec[z] += tmp2 - nu_as[z];
                                }
                                else{
                                    break;
                                }
                                
                            }
                        }

                        for (int z = 0; z < len_as; ++z)
                        {
                            tmp = tmpvec[z];
                            if(tmp>maxvals1[cord_spec(k, z, N)] ){
                                maxvals1[cord_spec(k, z, N)] = tmp;
                                if(ts[z] > twologn){
                                    maxvals2[cord_spec(k,z,N)] = tmp;
                                }
                            }
                        }

                        sort_k_largest_abs(cusum + cord_spec(0,v,p) , twologn, 0, p);
                        //partial_quicksort(cusum + cord_spec(0,j,p) , p , twologn);
                        double cumsum = 0;
                        int prev = 0;

                        for (int z = (len_as-1); z >= 0; --z)
                        {
                            if(ts[z] > twologn){
                                break;
                            }
                 
                            if(z==(len_as-1)){
                                prev = -1;
                            }
                            for (int hh = prev+1; hh < ts[z]; ++hh)
                            {
                                cumsum += cusum[cord_spec(hh,v,p)]*cusum[cord_spec(hh,v,p)];
                            }
                            prev = ts[z]-1;
                            tmp = cumsum;

                            if( (tmp ) >maxvals2[cord_spec(k,z,N)]){
                                maxvals2[cord_spec(k,z,N)] = tmp;
                            }
                        }
                    }
                }
            }
        }
    }

                
            
    
    SEXP thresholds1SEXP = PROTECT(allocVector(REALSXP, len_as ));
    double * thresholds1 = REAL(thresholds1SEXP);
    SEXP thresholds2SEXP = PROTECT(allocVector(REALSXP, len_as ));
    double * thresholds2 = REAL(thresholds2SEXP);

    for (int z = 0; z < len_as; ++z){
        sort_k_largest(maxvals1 + cord_spec(0,z,N) , toln, 0, N);
        sort_k_largest(maxvals2 + cord_spec(0,z,N) , toln, 0, N);
        thresholds1[z] = maxvals1[cord_spec(toln-1, z, N)];
        thresholds2[z] = maxvals2[cord_spec(toln-1, z, N)];

    }

    SEXP ret = PROTECT(allocVector(VECSXP, 4)); //the list to be returned in R
    SET_VECTOR_ELT(ret, 0, thresholds1SEXP);
    SET_VECTOR_ELT(ret, 1, thresholds2SEXP);
    SET_VECTOR_ELT(ret, 2, asI);
    SET_VECTOR_ELT(ret, 3, nu_asI);

    // creating list of names/titles to be returned in the output list
    SEXP names = PROTECT(allocVector(STRSXP, 4));
    //SET_STRING_ELT(names, 0, mkChar("CUSUM"));
    SET_STRING_ELT(names, 0, mkChar("without_partial"));
    SET_STRING_ELT(names, 1, mkChar("with_partial"));
    SET_STRING_ELT(names, 2, mkChar("as"));
    SET_STRING_ELT(names, 3, mkChar("nu_as"));


    setAttrib(ret, R_NamesSymbol, names);


    UNPROTECT(27);
    return(ret);











    







   /* SEXP out = PROTECT(allocVector(INTSXP, n));
    int * changepoints = INTEGER(out); //pointer to array
    //memset(changepoints, -1, sizeof(int)*n);
    for (size_t i = 0; i < n; i++) {
      changepoints[i] = -1;
    }
    int changepoint_counter = 0;
    int* changepoint_counter_ptr = &changepoint_counter;
    SEXP maxvalSEXP = PROTECT(allocVector(REALSXP, n));
    double * maxval = REAL(maxvalSEXP); //pointer to array
    memset(maxval, 0, sizeof(double)*n);
    SEXP startpointsSEXP = PROTECT(allocVector(INTSXP, n));
    int * startpoints = INTEGER(startpointsSEXP); //pointer to array
    memset(startpoints, 0, sizeof(int)*n);
    SEXP endpointsSEXP = PROTECT(allocVector(INTSXP, n));
    int * endpoints = INTEGER(endpointsSEXP); //pointer to array
    memset(endpoints, 0, sizeof(int)*n);
    SEXP maxaposesSEXP = PROTECT(allocVector(INTSXP, n));
    int * maxaposes = INTEGER(maxaposesSEXP); //pointer to array
    memset(maxaposes, 0, sizeof(int)*n);
    SEXP depthcounterSEXP= PROTECT(allocVector(INTSXP, n));
    int * depthcounter = INTEGER(depthcounterSEXP); //pointer to array
    memset(depthcounter, 0, sizeof(int)*n);
    SEXP coordschgSEXP = PROTECT(allocVector(INTSXP, n*p));
    int * coordchg = INTEGER(coordschgSEXP);
    memset(coordchg, 0, sizeof(int)*n*p);
    // first we compute all the cumulative sums of all
    // coordinates:

    




    // record all segment starts
    // we don't compute all cusums here but do it on the fly so we dont have to do anything
    // unecessary


    SEXP maxvaluesSEXP = PROTECT(allocVector(REALSXP, n*lenLens));
    double * maxvalues = REAL(maxvaluesSEXP);
    //memset(maxvalues, 0, sizeof(double)*n*lenLens);
    for (int i = 0; i < n*lenLens; ++i)
    {
        maxvalues[i] = NEGINF;
    }
    
    SEXP maxposSEXP = PROTECT(allocVector(INTSXP, n*lenLens));
    int * maxpos = INTEGER(maxposSEXP);
    memset(maxpos, 0, sizeof(int)*n*lenLens);
    SEXP maxasSEXP = PROTECT(allocVector(INTSXP, n*lenLens));
    int * maxas = INTEGER(maxasSEXP);
    memset(maxas, 0, sizeof(int)*n*lenLens);
    

    if(debug){
      for (int j = 0; j < lenLens; ++j)
      {
          //counter=0;



          for (int i = 0; i < n; i++)
          {
              //cord_spec(r,c, D) ((r) + (D)*(c))
              //compute_cusum(cumsum, i, i+len, &(maxpos[cord_spec(i,j,n)]), &(maxvalues[cord_spec(i,j,n)]));
              //segstarts[cord_spec(counter++,j,n)] = i;
              if(debug){
                  Rprintf("segstarts[%d, %d] = %d\n",i, j, segstarts[cord_spec(i,j,n)] );
              }


          }
      }
    }*/



}



SEXP cESAC_single(SEXP XI,SEXP nI, SEXP pI,SEXP thresholdsI,
    SEXP asI, SEXP nu_asI, SEXP len_asI, SEXP debugI){
    // X : p \times n
    PROTECT(XI);
    PROTECT(thresholdsI);
    
    PROTECT(asI);
    PROTECT(nu_asI);
    PROTECT(nI);
    PROTECT(pI);
    PROTECT(len_asI);
    PROTECT(debugI);

    double * X = REAL(XI);
    int n = *(INTEGER(nI));
    int p = *(INTEGER(pI));
    double *thresholds = (REAL(thresholdsI));
    
    double * as = REAL(asI);
    double * nu_as = REAL(nu_asI);
    int len_as = *(INTEGER(len_asI));
    int debug = *INTEGER(debugI);

    //UNPROTECT(3); // unprotecting all except the arrays
    if(debug){
      Rprintf("p = %d\n", p);
      Rprintf("n = %d\n", n);
    }


    SEXP tmpvecSEXP = PROTECT(allocVector(REALSXP, len_as));
    double * tmpvec = REAL(tmpvecSEXP);
    memset(tmpvec, 0, sizeof(double)*len_as);

    SEXP cumsumsSEXP = PROTECT(allocVector(REALSXP, p * (n+1)));
    double *cumsums = REAL(cumsumsSEXP); // p \times (n+1). first col is 0
    memset(cumsums, 0, p*(n+1)*sizeof(double));

    for (int j = 1; j <=n; ++j)
    {
        for (int i = 0; i < p; ++i)
        {
            //#define cord_spec(r,c, D) ((r) + (D)*(c))
            cumsums[cord_spec(i,j,p)] = X[cord_spec(i,j-1, p)] +cumsums[cord_spec(i,j-1, p)];
        }
    }

    SEXP cusumSEXP = PROTECT(allocVector(REALSXP, p*(n)));
    double * cusum = REAL(cusumSEXP);
    memset(cusum, 0, sizeof(double)*p*(n));
    int maxlen = p;
    int minlen = n;
    if(n>p){
        maxlen = n;
        minlen = p;
    }


    SEXP maxposSEXP = PROTECT(allocVector(INTSXP, 1));
    SEXP maxaposSEXP = PROTECT(allocVector(INTSXP, 1));
    int * maxpos = INTEGER(maxposSEXP);
    int * maxapos = INTEGER(maxaposSEXP);
    *maxpos = -10;
    double maximum = NEGINF;
    *maxapos = 0;
    int s = -1;
    int e = n-1;
    internal_find_changepoint( cumsums, cusum, maxpos,  &maximum, maxapos,
                             s, e, p, NULL, thresholds,
                            as, nu_as, len_as,  tmpvec,0, NULL,  debug);


    SEXP ret = PROTECT(allocVector(VECSXP, 2)); //the list to be returned in R
    SET_VECTOR_ELT(ret, 0, maxposSEXP);
    SET_VECTOR_ELT(ret, 1, maxaposSEXP);

    // creating list of names/titles to be returned in the output list
    SEXP names = PROTECT(allocVector(STRSXP, 2));
    //SET_STRING_ELT(names, 0, mkChar("CUSUM"));
    SET_STRING_ELT(names, 0, mkChar("pos"));
    SET_STRING_ELT(names, 1, mkChar("apos"));


    setAttrib(ret, R_NamesSymbol, names);

    UNPROTECT(15);
    return ret;
}


SEXP cESAC_test(SEXP XI,SEXP nI, SEXP pI,SEXP thresholdsI,
    SEXP asI, SEXP nu_asI, SEXP len_asI,  SEXP twolognI,SEXP tsI, 
    SEXP rescale_variance_boolI, SEXP fastI,SEXP debugI){

    PROTECT(XI);
    PROTECT(thresholdsI);
    
    PROTECT(asI);
    PROTECT(nu_asI);
    PROTECT(nI);
    PROTECT(pI);
    PROTECT(len_asI);
    PROTECT(twolognI);
    PROTECT(debugI);
    PROTECT(tsI);
    PROTECT(rescale_variance_boolI);
    PROTECT(fastI);
    double * X = REAL(XI);
    int n = *(INTEGER(nI));
    int p = *(INTEGER(pI));
    double *thresholds = (REAL(thresholdsI));
    
    double * as = REAL(asI);
    double * nu_as = REAL(nu_asI);
    int len_as = *(INTEGER(len_asI));
    int twologn = *INTEGER(twolognI);
    int *ts = INTEGER(tsI);
    int debug = *INTEGER(debugI);
    int fast = *INTEGER(fastI);
    int rescale_variance_bool = *INTEGER(rescale_variance_boolI);

    int maxlen = p;
    int minlen = n;
    if(n>p){
        maxlen = n;
        minlen = p;
    }
    SEXP vectorSEXP = PROTECT(allocVector(REALSXP, maxlen));

    double * vector = REAL(vectorSEXP);
    memset(vector, 0, sizeof(double)*maxlen);

    
    if (rescale_variance_bool)
    {
        rescale_variance(X, NULL, n, p, vector);
    }

    SEXP cumsumsSEXP = PROTECT(allocVector(REALSXP, p * (n+1)));
    double *cumsums = REAL(cumsumsSEXP); // p \times (n+1). first col is 0
    memset(cumsums, 0, p*(n+1)*sizeof(double));

    for (int j = 1; j <=n; ++j)
    {
        for (int i = 0; i < p; ++i)
        {
            //#define cord_spec(r,c, D) ((r) + (D)*(c))
            cumsums[cord_spec(i,j,p)] = X[cord_spec(i,j-1, p)] +cumsums[cord_spec(i,j-1, p)];
        }
    }

    int lencusum = p*n;
    if(fast){lencusum = p;}
    SEXP cusumSEXP = PROTECT(allocVector(REALSXP, lencusum));
    double * cusum = REAL(cusumSEXP);
    memset(cusum, 0, sizeof(double)*lencusum);

    double tmp = 0;
    double tmpmax = 0;
    int tmppos = 0;
    int tmpapos = 0;
    //double a_tmp=-100000;
    //int pos_a_tmp = 0;
    //double * val=0;
    //int n = e-s;

    //int detected = 0;
    double tmp2 =0;
    int localdetected = 0;
    int end = 1;

    SEXP detectedSEXP = PROTECT(allocVector(INTSXP,1));
    int * detected = INTEGER(detectedSEXP);
    *detected = 0;

    SEXP tmpvecSEXP = PROTECT(allocVector(REALSXP, len_as));
    double * tmpvec = REAL(tmpvecSEXP);
    memset(tmpvec, 0, sizeof(double)*len_as);

    //int detected = 0;
    int j = 0;
    int s = -1;
    int e = n-1;
    int k = 1;
    if(fast){
        k = (s+e)/2;

        singleCUSUM(cumsums, cusum, s, e, p, k);
        memset(tmpvec, 0, sizeof(double)*len_as);
        
        // first aggregate thresholded CUSUMs
        j = 0;
        for (int i = 0; i < p; ++i){
            tmp2 = cusum[cord_spec(i,j, p)]*cusum[cord_spec(i,j, p)];
            for (int z = 0; z < len_as; ++z){
                if (fabs(cusum[cord_spec(i,j, p)])>as[z])
                {
                    tmpvec[z] += tmp2 - nu_as[z];
                }
                else{
                    break;
                }
                
            }
        }

        // check if T >0
        for (int z = 0; z < len_as; ++z)
        {
            tmp = tmpvec[z] - thresholds[z];
            if(ts[z] <= twologn){
                break;
            }
            else if(tmp>0){
                *detected = 1;
                break;
            }
        }

        // if T<= 0, check if P >0:
        if(twologn>0 && *detected ==0){
            sort_k_largest_abs(cusum + cord_spec(0,j,p) , twologn, 0, p);
            double cumsum = 0;
            int prev = 0;

            for (int z = (len_as-1); z >= 0; --z)
            {
                if(ts[z] > twologn){
                    break;
                }
     
                if(z==(len_as-1)){
                    prev = -1;
                }
                for (int hh = prev+1; hh < ts[z]; ++hh)
                {
                    cumsum += cusum[cord_spec(hh,j,p)]*cusum[cord_spec(hh,j,p)];
                }
                prev = ts[z]-1;
                tmp = cumsum;

                if( (tmp ) >thresholds[z]){
                    *detected = 1;
                    break;
                }
            }
        }

    }
    else{

        //CUSUM(cumsums, cusum, s, e, p);
        end = e-s-1;
        for (j = 0; j < end; ++j){
            singleCUSUM(cumsums, cusum, s, e, p, j);
            memset(tmpvec, 0, sizeof(double)*len_as);
            
            // first aggregate thresholded CUSUMs
            for (int i = 0; i < p; ++i){
                tmp2 = cusum[cord_spec(i,0, p)]*cusum[cord_spec(i,0, p)];
                for (int z = 0; z < len_as; ++z){
                    if (fabs(cusum[cord_spec(i,0, p)])>as[z])
                    {
                        tmpvec[z] += tmp2 - nu_as[z];
                    }
                    else{
                        break;
                    }
                    
                }
            }



            // check if T >0

            if(!*detected){
                for (int z = 0; z < len_as; ++z)
                {
                    /*tmp = tmpvec[z] - thresholds[z];
                    if (tmp> *maximum)
                    {
                        *maximum= tmp;
                        *maxpos = s+j+1;
                        *maxa_pos = z;
                    }*/
                    tmp = tmpvec[z] - thresholds[z];
                    if(tmp>0 && ts[z] > twologn){
                        *detected = 1;
                        break;
                    }
                }
                if(*detected){
                    break;
                }

                // if T<= 0, check if P >0:
                if(twologn>0){
                    sort_k_largest_abs(cusum + cord_spec(0,0,p) , twologn, 0, p);
                    //partial_quicksort(cusum + cord_spec(0,j,p) , p , twologn);
                    double cumsum = 0;
                    int prev = 0;

                    for (int z = (len_as-1); z >= 0; --z)
                    {
                        if(ts[z] > twologn){
                            break;
                        }
             
                        if(z==(len_as-1)){
                            prev = -1;
                        }
                        for (int hh = prev+1; hh < ts[z]; ++hh)
                        {
                            cumsum += cusum[cord_spec(hh,0,p)]*cusum[cord_spec(hh,0,p)];
                        }
                        prev = ts[z]-1;
                        tmp = cumsum;

                        if(tmp >thresholds[z]){
                            *detected = 1;
                            break;
                        }
                    }
                }
            
                if(*detected){
                    break;
                }
            }
        }

    }




    UNPROTECT(17);
    return(detectedSEXP);


}






SEXP cESAC_test_calibrate(SEXP nI, SEXP pI, SEXP NI, SEXP tolnI, SEXP asI, SEXP nu_asI, SEXP len_asI, 
                          SEXP twolognI,SEXP tsI, SEXP fastI, SEXP rescale_variance_boolI, SEXP debugI){

    PROTECT(nI);
    PROTECT(pI);
    PROTECT(NI);
    PROTECT(tolnI);
    PROTECT(asI);
    PROTECT(nu_asI);
    PROTECT(len_asI);
    PROTECT(twolognI);
    PROTECT(tsI);
    PROTECT(debugI);
    PROTECT(fastI);
    PROTECT(rescale_variance_boolI);
    int n = *(INTEGER(nI));
    int p = *(INTEGER(pI));
    int N = *(INTEGER(NI));
    int toln = *(INTEGER(tolnI));
    double * as = REAL(asI);
    double * nu_as = REAL(nu_asI);
    int len_as = *(INTEGER(len_asI));
    int twologn = *INTEGER(twolognI);
    int *ts = INTEGER(tsI);
    int debug = *INTEGER(debugI);
    int fast = *INTEGER(fastI);
    int rescale_variance_bool = *INTEGER(rescale_variance_boolI);

    if(debug){
        Rprintf("twologn = %d\n", twologn);
    }

    // SEXP Xsexp = PROTECT(allocVector(REALSXP, n*p));
    // double * X = REAL(Xsexp);

    // max vals for no partial sum
    SEXP maxvals1SEXP = PROTECT(allocVector(REALSXP, len_as * N));
    double * maxvals1 = REAL(maxvals1SEXP);
    //memset(maxvals1, 0, sizeof(double)*len_as*N);

    // max vals for partial sum
    SEXP maxvals2SEXP = PROTECT(allocVector(REALSXP, len_as * N));
    double * maxvals2 = REAL(maxvals2SEXP);

    int maxlen = p;
    int minlen = n;
    if(n>p){
        maxlen = n;
        minlen = p;
    }
    SEXP vectorSEXP = PROTECT(allocVector(REALSXP, maxlen));

    double * vector = REAL(vectorSEXP);
    memset(vector, 0, sizeof(double)*maxlen);

    for (int i = 0; i < len_as * N; ++i)
    {
        maxvals2[i] = NEGINF;
        maxvals1[i] = NEGINF;
    }
    //memset(maxvals2, 0, sizeof(double)*len_as*N);

    SEXP testposSEXP = PROTECT(allocVector(INTSXP, n));
    int *testpos  = INTEGER(testposSEXP);
    memset(testpos, 0, sizeof(int)*n);

    // find test poses
    int count = 0;
    int s = -1;
    int e = n-1;
    int middle = (s+e)/2;

    if(fast){
        count = 1;
        testpos[0] = middle;
    }
    else{
        count = n-1;
        for (int i = 0; i < count; ++i)
        {
            testpos[i] = i;
        }
    }

/*    testpos[count++] = middle;
    int jump = 2;

    while(middle + jump < e || middle - jump > s){
        if(middle + jump <e){
            testpos[count++] = middle + jump;
        }
        if(middle - jump > s){
            testpos[count++] = middle - jump;
        }
        jump = jump*2;
    }

    R_qsort_int(testpos,1,count);*/

    SEXP tmpvecSEXP = PROTECT(allocVector(REALSXP, len_as));
    double * tmpvec = REAL(tmpvecSEXP);
    //memset(tmpvec, 0, sizeof(double)*len_as);
    
    int testpoint;
    double tmp2;
    double tmp;
  
    SEXP XSEXP = PROTECT(allocVector(REALSXP, p * (n)));
    double *X = REAL(XSEXP);
    memset(X, 0, sizeof(double)*p*n);
  
    SEXP cumsumsSEXP = PROTECT(allocVector(REALSXP, p * (n+1)));
    double *cumsums = REAL(cumsumsSEXP); // p \times (n+1). first col is 0
    memset(cumsums, 0, sizeof(double)*p*(n+1));
    
    SEXP cusumSEXP = PROTECT(allocVector(REALSXP, p));
    double * cusum = REAL(cusumSEXP);
    memset(cusum, 0, sizeof(double)*p);
    
    int i=0;
    int j = 0;

    for (int k = 0; k < N; ++k)
    {
        
        // generate X
        GetRNGstate();
        for (j = 0; j < n; ++j)
        {
            for (i = 0; i < p; ++i)
            {
                X[cord_spec(i,j,p)] = norm_rand();
            }
        }
        PutRNGstate();

        // normalize / rescale variance
        if (rescale_variance_bool)
        {
            rescale_variance(X,NULL, n, p, vector);
        }
        
        memset(cumsums, 0, p*sizeof(double));

        for (i = 0; i < p; ++i)
        {
            for (j = 0; j < n; ++j)
            {
                cumsums[cord_spec(i,j+1,p)] = X[cord_spec(i,j,p)] +cumsums[cord_spec(i,j, p)];
            }
        }

        for (j = 0; j < count; ++j)
        {
            testpoint = testpos[j];
            singleCUSUM(cumsums, cusum, s, e, p, testpoint);
            memset(tmpvec, 0, sizeof(double)*len_as);
            
            // first aggregate thresholded CUSUMs
            for (i = 0; i < p; ++i){
                tmp2 = cusum[cord_spec(i,0, p)]*cusum[cord_spec(i,0, p)];
                for (int z = 0; z < len_as; ++z){
                    if (fabs(cusum[cord_spec(i,0, p)])>as[z])
                    {
                        tmpvec[z] += tmp2 - nu_as[z];
                    }
                    else{
                        break;
                    }
                    
                }
            }

            for (int z = 0; z < len_as; ++z)
            {
                tmp = tmpvec[z];
                if(tmp > maxvals1[cord_spec(k,z, N)]){
                    maxvals1[cord_spec(k,z, N)] = tmp;
                    if(ts[z] > twologn){
                        maxvals2[cord_spec(k,z, N)] = tmp;
                    }
                }
            }

            if(twologn>0){
                sort_k_largest_abs(cusum , twologn, 0, p);
                //partial_quicksort(cusum + cord_spec(0,j,p) , p , twologn);
                double cumsum = 0;
                int prev = 0;

                for (int z = (len_as-1); z >= 0; --z)
                {
                    if(ts[z] > twologn){
                        break;
                    }
         
                    if(z==(len_as-1)){
                        prev = -1;
                    }
                    for (int hh = prev+1; hh < ts[z]; ++hh)
                    {
                        cumsum += cusum[cord_spec(hh,0,p)]*cusum[cord_spec(hh,0,p)];
                    }
                    prev = ts[z]-1;
                    tmp = cumsum ;
                    //Rprintf("%d max = %f\n", ts[z], tmp);

                    if(tmp>maxvals2[cord_spec(k,z, N)]){
                        maxvals2[cord_spec(k, z,N)] = tmp;
                    }
                }
            }
            
          

        }
        




    }

    
    SEXP thresholds1SEXP = PROTECT(allocVector(REALSXP, len_as ));
    double * thresholds1 = REAL(thresholds1SEXP);
    SEXP thresholds2SEXP = PROTECT(allocVector(REALSXP, len_as ));
    double * thresholds2 = REAL(thresholds2SEXP);

    for (int z = 0; z < len_as; ++z){
        //sort_k_largest(maxvals1 + cord_spec(0,z,N) , toln, 0, N);
        R_qsort(maxvals1 + cord_spec(0,z,N) , 1,N);
        //sort_k_largest(maxvals2 + cord_spec(0,z,N) , toln, 0, N);
        R_qsort(maxvals2 + cord_spec(0,z,N) , 1,N);
        //thresholds1[z] = maxvals1[cord_spec(toln-1, z, N)];
        thresholds1[z] = maxvals1[cord_spec(N-toln, z, N)];
        //thresholds2[z] = maxvals2[cord_spec(toln-1, z, N)];
        thresholds2[z] = maxvals2[cord_spec(N-toln, z, N)];

    }

    SEXP ret = PROTECT(allocVector(VECSXP, 4)); //the list to be returned in R
    SET_VECTOR_ELT(ret, 0, thresholds1SEXP);
    SET_VECTOR_ELT(ret, 1, thresholds2SEXP);
    SET_VECTOR_ELT(ret, 2, asI);
    SET_VECTOR_ELT(ret, 3, nu_asI);

    // creating list of names/titles to be returned in the output list
    SEXP names = PROTECT(allocVector(STRSXP, 4));
    //SET_STRING_ELT(names, 0, mkChar("CUSUM"));
    SET_STRING_ELT(names, 0, mkChar("without_partial"));
    SET_STRING_ELT(names, 1, mkChar("with_partial"));
    SET_STRING_ELT(names, 2, mkChar("as"));
    SET_STRING_ELT(names, 3, mkChar("nu_as"));


    setAttrib(ret, R_NamesSymbol, names);


    UNPROTECT(24);
    return(ret);


}

