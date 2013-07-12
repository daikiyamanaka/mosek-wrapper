#include "QuadricOptimizer.h"
#include "mosek.h"
QuadricOptimizer::QuadricOptimizer(){

}

QuadricOptimizer::~QuadricOptimizer(){
	
}

static void MSKAPI printstr(void *handle,
                            MSKCONST char str[])
{
  printf("%s",str);
} /* printstr */

void QuadricOptimizer::solve(Eigen::SparseMatrix<double> &Q, Eigen::VectorXd &C, Eigen::SparseMatrix<double> &A, Eigen::VectorXd &b, Eigen::VectorXd &x){
    bool optimized = false;
    NUMCON = A.rows();
    NUMVAR = Q.rows();
    NUMANZ = A.nonZeros();

	int num_qnz = 0;
	for (int k=0; k<Q.outerSize(); ++k){		     	
		for (Eigen::SparseMatrix<double>::InnerIterator it(Q,k); it; ++it)
		{
			if(it.col() <= it.row()){					
				num_qnz ++;
				if(it.row() < 0 || it.col() < 0){
					std::cout << "index is minus " << it.row() << " " << it.col() << std::endl;
				}
			}
		}
	}
	NUMQNZ = num_qnz;

  std::cout << "NUMCON: " << NUMCON << std::endl;
  std::cout << "NUMVAR: " << NUMVAR << std::endl;
  std::cout << "NUMANZ: " << NUMANZ << std::endl;
  std::cout << "NUMQNZ: " << NUMQNZ << std::endl;  

  MSKrescodee  r;

  double       *c    = new double[NUMVAR];

  /*
  MSKboundkeye bkc[]  = {MSK_BK_FX, MSK_BK_FX};
  double       blc[]  = {16.0/(double)width, 32.0/(double)height};
  double       buc[]  = {16.0/(double)width, 32.0/(double)height}; 
  */
  // MSKboundkeye bkc[]  = {MSK_BK_FX, MSK_BK_FX, MSK_BK_FX};
  // double       blc[]  = {0.0, 0.0, 1.0};
  // double       buc[]  = {0.0, 0.0, 1.0}; 
  MSKboundkeye *bkc = new MSKboundkeye[b.size()];
  double *blc = new double[b.size()];
  double *buc = new double[b.size()];  

  for(int i=0; i<b.size(); i++){
  	bkc[i] = MSK_BK_FX;
  	blc[i] = b(i);
  	buc[i] = b(i);
  }  

  MSKboundkeye *bkx  = new MSKboundkeye[NUMVAR];	  
  double *blx = new double[NUMVAR];
  double *bux = new double[NUMVAR];
  for(int i=0; i<NUMVAR; i++){
    bkx[i] = MSK_BK_RA;
    //blx[i] = C(i)/(-2.0);
    blx[i] = 0.0;
  	bux[i] = 1.0;
    c[i] = C(i);
  }

  // condition matrix 
  MSKint32t    *aptrb = new MSKint32t[A.cols()];//number of column
  MSKint32t    *aptre = new MSKint32t[A.cols()]; //number of column    
  MSKint32t    *asub = new MSKint32t[A.nonZeros()];// row index
  double       *aval = new double[A.nonZeros()];

  //std::cout << A.nonZeros() << std::endl;
  int index = 0;
  for (int k=0; k<A.outerSize(); ++k){
  	aptrb[k] = index;
  	for (Eigen::SparseMatrix<double>::InnerIterator it(A,k); it; ++it)
  	{
  		aval[index] = it.value();    		
  		asub[index] = it.row();
      index ++;
    }
    aptre[k] = index;
  }


  //square matrix    
  MSKint32t    *qsubi = new MSKint32t[NUMQNZ];
  MSKint32t    *qsubj = new MSKint32t[NUMQNZ];    
  double       *qval = new double[NUMQNZ];

  MSKint32t    j,i;
  double       *xx = new double[NUMVAR];
  MSKenv_t     env;
  MSKtask_t    task;

    /* Create the mosek environment. */
  r = MSK_makeenv(&env,NULL);
  if ( r==MSK_RES_OK )
  {   
  /* Create the optimization task. */
	r = MSK_maketask(env,NUMCON,NUMVAR,&task);
  	
	if ( r==MSK_RES_OK )
	{
		r = MSK_linkfunctotaskstream(task,MSK_STREAM_LOG,NULL,printstr);

  /* Append 'NUMCON' empty constraints. The constraints will initially have no bounds. */
		if ( r == MSK_RES_OK )
			r = MSK_appendcons(task,NUMCON); 

      /* Append 'NUMVAR' variables. The variables will initially be fixed at zero (x=0). */
		if ( r == MSK_RES_OK )
			r = MSK_appendvars(task,NUMVAR);

      /* Optionally add a constant term to the objective. */
		if ( r ==MSK_RES_OK )
			r = MSK_putcfix(task,0.0);

		for(j=0; j<NUMVAR && r == MSK_RES_OK; ++j)
		{
	  /* Set the linear term c_j in the objective.*/  
			if(r == MSK_RES_OK)
				r = MSK_putcj(task,j,c[j]);
    			
	  /* Set the bounds on variable j. blx[j] <= x_j <= bux[j] */
			if(r == MSK_RES_OK)
				r = MSK_putvarbound(task,
				j,           /* Index of variable.*/
				bkx[j],      /* Bound key.*/
				blx[j],      /* Numerical value of lower bound.*/
				bux[j]);     /* Numerical value of upper bound.*/
    					
	  /* Input column j of A */   
			if(r == MSK_RES_OK)
				r = MSK_putacol(task,
			    j,                 /* Variable (column) index.*/
			    aptre[j]-aptrb[j], /* Number of non-zeros in column j.*/
			    asub+aptrb[j],     /* Pointer to row indexes of column j.*/
			    aval+aptrb[j]);    /* Pointer to Values of column j.*/
    
		}

      /* Set the bounds on constraints. for i=1, ...,NUMCON : blc[i] <= constraint i <= buc[i] */
		for(i=0; i<NUMCON && r==MSK_RES_OK; ++i)
			r = MSK_putconbound(task,
			    i,           /* Index of constraint.*/
			    bkc[i],      /* Bound key.*/
			    blc[i],      /* Numerical value of lower bound.*/
			    buc[i]);     /* Numerical value of upper bound.*/

			if ( r==MSK_RES_OK )
			{
	    /*
	     * The lower triangular part of the Q^o
	     * matrix in the objective is specified.
	     */

	      // qsubi[0] = 0;   qsubj[0] = 0;  qval[0] = 2.0;
	      // qsubi[1] = 1;   qsubj[1] = 1;  qval[1] = 0.2;
	      // qsubi[2] = 2;   qsubj[2] = 0;  qval[2] = -1.0;
	      // qsubi[3] = 2;   qsubj[3] = 2;  qval[3] = 2.0;
	     int obj_index = 0;
	     for (int k=0; k<Q.outerSize(); ++k){		     	
	     	 for (Eigen::SparseMatrix<double>::InnerIterator it(Q,k); it; ++it)
	     	 {
		     	 	 if(it.col() <= it.row()){
		     	 	 	qsubi[obj_index] = it.row();
		     	 	 	qsubj[obj_index] = it.col();
		     	 	 	qval[obj_index] = it.value();
		     	 	 	obj_index ++;
		     	 	 	if(it.row() < 0 || it.col() < 0){
		     	 	 		std::cout << "index is minus " << it.row() << " " << it.col() << std::endl;
		     	 	 	}
		     	 	 }
		     	 	}
	     	 }

	  /* Input the Q^o for the objective. */

	      r = MSK_putqobj(task,NUMQNZ,qsubi,qsubj,qval);
	    }

	    if ( r==MSK_RES_OK )
	    {
	  /*
	   * The lower triangular part of the Q^0
	   * matrix in the first constraint is specified.
       This corresponds to adding the term
       - x_1^2 - x_2^2 - 0.1 x_3^2 + 0.2 x_1 x_3
	  */        
        //qsubi[0] = 0;   qsubj[0] = 0;  qval[0] = -2.0;
        //qsubi[1] = 1;   qsubj[1] = 1;  qval[1] = -2.0;
        //qsubi[2] = 2;   qsubj[2] = 2;  qval[2] = -0.2;
        //qsubi[3] = 2;   qsubj[3] = 0;  qval[3] = 0.2;

	  /* Put Q^0 in constraint with index 0. */
        //r = MSK_putqconk(task, 0, 4, qsubi, qsubj, qval); 
      }

      if ( r==MSK_RES_OK )
      	r = MSK_putobjsense(task, MSK_OBJECTIVE_SENSE_MINIMIZE);      

      if ( r==MSK_RES_OK )
      {
      	MSKrescodee trmcode;
      
	  /* Run optimizer */
      	r = MSK_optimizetrm(task,&trmcode);

	  /* Print a summary containing information
	     about the solution for debugging purposes*/
      	MSK_solutionsummary (task,MSK_STREAM_LOG);
      
      	if ( r==MSK_RES_OK )
      	{
      		MSKsolstae solsta;
      		int j;

      		MSK_getsolsta (task,MSK_SOL_ITR,&solsta);

      		switch(solsta)
      		{
      			case MSK_SOL_STA_OPTIMAL:   
      			case MSK_SOL_STA_NEAR_OPTIMAL:
      			MSK_getxx(task,
			        MSK_SOL_ITR,    /* Request the interior solution. */
      				xx);
            
      			printf("Optimal primal solution\n");
      			for(j=0; j<NUMVAR; ++j)
      				//printf("x[%d]: %e\n",j,xx[j]);
                optimized = true;
      			break;
      			case MSK_SOL_STA_DUAL_INFEAS_CER:
      			case MSK_SOL_STA_PRIM_INFEAS_CER:
      			case MSK_SOL_STA_NEAR_DUAL_INFEAS_CER:
      			case MSK_SOL_STA_NEAR_PRIM_INFEAS_CER:  
      			printf("Primal or dual infeasibility certificate found.\n");
      			break;

      			case MSK_SOL_STA_UNKNOWN:
      			printf("The status of the solution could not be determined.\n");
      			break;
      			default:
      			printf("Other solution status.");
      			break;
      		}
      	}
      	else
      	{
      		printf("Error while optimizing.\n");
      	}
      }

      if (r != MSK_RES_OK)
      {
	  /* In case of an error print error code and description. */      
      	char symname[MSK_MAX_STR_LEN];
      	char desc[MSK_MAX_STR_LEN];

      	printf("An error occurred while optimizing.\n");     
      	MSK_getcodedesc (r,
      		symname,
      		desc);
      	printf("Error %s - '%s'\n",symname,desc);
      }
    }       
    MSK_deletetask(&task);
  }
  MSK_deleteenv(&env);
  if(optimized){
    x = Eigen::VectorXd(NUMVAR);
    for(int i=0; i<NUMVAR; i++){
        x(i) = xx[i];
    }
  }
}

