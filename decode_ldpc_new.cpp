
/*LDPC Decoder*/

#include "mex.h"
#include "matrix.h"			// for Matlab mx and mex fuctions
#include "math.h"
#include <stdlib.h>			// what for
#include "decodeutil_new.h"

#define BIGVALUE_COLS -1
#define BIGVALUE_ROWS -1


// tried out the precision limit and this are good reluts in term of 
// precision and velocity
double lntanh(double x)
{
    if(x > 37) return 1.1102*pow(10,-16);
    
    else if(x < 1.1102*pow(10,-16)) return 37.43;
    
    else return -log(tanh(x/2));    
}

void decode(double max_iter, double *vhat, int mrows, int ncols, double *iter_out, double *LLR,
        double *check_node_ones,double max_check_degree,
        double *variable_node_ones,double max_variable_degree)
        
{
    double **variable_messages,**check_messages, **bitmessage_temp ;
    double **check_node_ones_matrix,**variable_node_ones_matrix;
    
    variable_messages = matrix(0,mrows-1,0,ncols-1);
    check_messages = matrix(0,mrows-1,0,ncols-1);
    
    check_node_ones_matrix = matrix(0,mrows-1,0,max_check_degree-1);
    variable_node_ones_matrix = matrix(0,max_variable_degree-1,0,ncols-1);
    
    
    
    for (int i = 0; i < max_check_degree; i++)
        for(int j = 0; j< mrows; j++)
            check_node_ones_matrix[j][i] =  *check_node_ones++;
    
    for(int i = 0; i < ncols; i++)
        for(int j = 0; j < max_variable_degree; j++)
            variable_node_ones_matrix[j][i] = *variable_node_ones++;
   
    double element = 0, col_index = 0, row_index = 0;
    
    //initializing the matrices
    for (int i = 0; i < ncols; i++)
    {
        for(int j = 0; j < max_variable_degree; j++)
        {
            row_index = variable_node_ones_matrix[j][i];
            if(row_index == BIGVALUE_ROWS)
            {
                break;
            }
            
            variable_messages[(int)row_index][i] = LLR[i];
        }
    }
    
    int iteration;
    double temp = 0, prodi = 1, temp_2, temp_3;
    int parity;
    double *sum_of_b;
    sum_of_b = vector(0,ncols-1);
    double ** temp_matrix;
    temp_matrix = matrix(0,mrows-1,0,ncols-1);
    
    for (iteration = 0; iteration < max_iter; iteration++)
    {   
        
        // For speed-up the program
        for(int u = 0; u < mrows; u++)
        {
            prodi = 1;
            temp = 0;
            for(int v = 0; v < max_check_degree; v++)
            {
                col_index = check_node_ones_matrix[u][v];
                if(col_index != BIGVALUE_COLS)
                {
                    temp_2 = variable_messages[u][(int)col_index];
                    if(temp_2 < 0)
                        prodi = prodi * (-1);

                    temp = temp + lntanh(abs(temp_2));
                }
            }
            
            for(int v = 0; v < max_check_degree; v++)
            {
                col_index = check_node_ones_matrix[u][v];
                if(col_index != BIGVALUE_COLS)
                    temp_matrix[u][(int) col_index] = prodi * temp;
            }
        }
        
        for (int u = 0; u < mrows; u++)
        {
            for (int v = 0; v < max_check_degree; v++)
            {
                temp_2 = 0;
                temp = 0;
                temp_3 = 0;
                prodi = 1;
                
                element = check_node_ones_matrix[u][v];
                if (element == BIGVALUE_COLS)
                    break;
                
                
                temp_3 = temp_matrix[u][(int)element];
                if(temp_3 < 0)
                {
                    prodi *= (-1);
                    temp_3 *= (-1);
                }
                
                temp_2 = variable_messages[u][(int)element];
                
                if(temp_2 < 0)
                    prodi *= (-1);
                
                temp = temp_3 - lntanh(abs(temp_2));
                
                check_messages[u][(int)element] = lntanh(abs(temp)) * prodi;
            }
        }
         
        //sum across columns and hard decision
        for (int i = 0; i < ncols; i++)
        {
            temp = 0;
            for(int j = 0; j < max_variable_degree; j++)
            {
                row_index = variable_node_ones_matrix[j][i];
                if(row_index == BIGVALUE_ROWS)
                    break;
                
                temp = temp + check_messages[(int)row_index][i];
            }
            
            sum_of_b[i] = temp;
            vhat[i] = (temp + LLR[i] < 0 ? 1.0 : 0.0); 
        }
        
        // check if is a valid codeword
        for (int j = 0; j < mrows; j++)
        {
            parity = 0;
            for(int c = 0; c < max_check_degree; c++)
            {
                col_index = check_node_ones_matrix[j][c];
                if(col_index == -1)
                    break;
                parity = parity + (int) vhat[(int)col_index];
            }
            
            parity = parity % 2;
            
            if(parity == 1)
                break;
        }
        
        if (parity == 0)
        {
            *iter_out = iteration + 1;
            return;
        }
        
        //update the APP LLR 
        for (int i = 0; i < ncols; i++)
        {
            for(int j = 0; j < max_variable_degree; j++)
            {
                temp = 0;
                element = variable_node_ones_matrix[j][i];
                if(element == BIGVALUE_ROWS)
                    break;
                
                temp = sum_of_b[i] - check_messages[(int)element][i];
                
                variable_messages[(int)element][i] = temp + LLR[i];
            }          
        }
    }
    *iter_out = iteration;
}



void mexFunction( int nlhs, mxArray *plhs[],
        int nrhs, const mxArray*prhs[] )
{
    double *vhat, *iter_out,*LLR;
    double *check_node_ones, *variable_node_ones; /*pointer variables for input Matrices*/
    double max_iter,max_check_degree,max_variable_degree;
    int mrows,ncols;
    
    
    max_iter = mxGetScalar(prhs[0]); // maximum iterations
    
    LLR  = mxGetPr(prhs[1]); //pointer to initial APP LLR
    check_node_ones = mxGetPr(prhs[2]);  //pointer to matrix containing column indeces in each row which are '1'
    max_check_degree = mxGetScalar(prhs[3]); 
    
    variable_node_ones = mxGetPr(prhs[4]); //pointer to matrix containing row indeces in each column which are '1'
    max_variable_degree = mxGetScalar(prhs[5]);
    
    mrows = mxGetScalar(prhs[6]); // number of rows of H
    ncols = mxGetScalar(prhs[7]); // number of cols of H
    
    plhs[0] = mxCreateDoubleMatrix(ncols, 1, mxREAL); /*matrix for output*/
    vhat = mxGetPr(plhs[0]);	/*pointer to output*/
    
    plhs[1] = mxCreateDoubleScalar(0);
    iter_out = mxGetPr(plhs[1]);
    
    decode(max_iter,vhat,mrows,ncols,iter_out,LLR,check_node_ones,max_check_degree,variable_node_ones,max_variable_degree);
}
