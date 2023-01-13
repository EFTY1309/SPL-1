#include<bits/stdc++.h>
using namespace std;
#define N 3

void transpose(long double W[N][N],long double transpose[N][N], int row,int col)
{
    for (int i = 0; i < row; i++)
    {
        for (int j = 0; j < col; j++)
        {
            transpose[j][i] = W[i][j];
        }
    }
}


void multiplyMatrices(long double firstMatrix[N][N],long  double secondMatrix[N][N],long  double mult[N][N], int rowFirst, int columnFirst, int rowSecond, int columnSecond)
{
	int i, j, k;
	for(i = 0; i < rowFirst; ++i)
	{
		for(j = 0; j < columnSecond; ++j)
		{
			mult[i][j] = 0;
		}
	}

	for(i = 0; i < rowFirst; ++i)
	{
		for(j = 0; j < columnSecond; ++j)
		{
			for(k=0; k<columnFirst; ++k)
			{
				mult[i][j] += firstMatrix[i][k] * secondMatrix[k][j];
			}
		}
	}
}

double frobeniusnorm(long double y[N][N],int row,int col)
{
    double sum=0,norm;

    for(int i=0;i<row;i++)
        for(int j=0;j<col;j++)
    {
        sum+=pow(y[i][j],2);
    }
    norm=sqrt(sum);
    return norm;
}
void display(long double a[N][N],int row,int col)
{
    for(int i=0;i<row;i++){
      for(int j=0;j<col;j++)
      cout<<a[i][j]<<"  ";
      cout<<endl;
    }
}

void ElementwiseMultiplication(long double a[N][N],long double b[N][N],int row,int col)
{
    for(int i=0;i<row;i++)
        for(int j=0;j<col;j++)
    {
        a[i][j]=(a[i][j])*(b[i][j]);
    }


}





int main(void)
{
 long double x[N][N];
  int row,col;
  cout<<"Enter the number of row and column for X matrix:"<<endl;
  cin>>row>>col;

  for(int i=0;i<row;i++)
    for(int j=0;j<col;j++)
  {
      cin>>x[i][j];
  }
  long double a[N][N];
  int row1,col1;

  cout<<"Enter the number of row and column for A matrix:"<<endl;
  cin>>row1>>col1;

  for(int i=0;i<row1;i++)
    for(int j=0;j<col1;j++)
  {
      cin>>a[i][j];
  }

 long  double b[N][N];
  int row2,col2;

  cout<<"Enter the number of row and column for B matrix:"<<endl;
  cin>>row2>>col2;

  for(int i=0;i<row2;i++)
    for(int j=0;j<col2;j++)
  {
      cin>>b[i][j];
  }

  int t=7;
  while(t--)
  {
  long  double a_t[N][N],b_t[N][N],xb_t[N][N],ab[N][N],abb_t[N][N],a_tx[N][N],a_ta[N][N] ,a_tab[N][N],ab1[N][N],L[N][N];
  long double a1[N][N],b1[N][N];
  transpose(a,a_t,row1,col1);
  transpose(b,b_t,row2,col2);

  multiplyMatrices(x,b_t,xb_t,row,col,col2,row2);
  multiplyMatrices(a,b,ab,row1,col1,row2,col2);
  multiplyMatrices(ab,b_t,abb_t,row1,col2,col2,row2);



  for(int i=0;i<row1;i++){
    for(int j=0;j<col1;j++){
        a1[i][j]=(xb_t[i][j])/(abb_t[i][j]);
    }
  }
  cout<<endl<<"Updated a1 matrix is:"<<endl;
  display(a1,row1,col1);
  cout<<endl;

  ElementwiseMultiplication(a,a1,row1,col1);

  multiplyMatrices(a_t,x,a_tx,col1,row1,row,col);
  multiplyMatrices(a_t,a,a_ta,col1,row1,row1,col1);
  multiplyMatrices(a_ta,b,a_tab,col1,col1,row2,col2);


  cout<<"updated A matrix is:"<<endl;
  display(a,row1,col1);
  cout<<endl;

  for(int i=0;i<row2;i++){
    for(int j=0;j<col2;j++){
        b1[i][j]=(a_tx[i][j])/(a_tab[i][j]);
    }
    cout<<endl;
  }
  cout<<endl<<"Updated b1 matrix is:"<<endl;
  display(b1,row2,col2);
  cout<<endl;


  ElementwiseMultiplication(b,b1,row2,col2);

  cout<<"updated B matrix is:"<<endl;

   display(b,row2,col2);

  multiplyMatrices(a,b,ab1,row1,col1,row2,col2);
  cout<<"Iteration:"<<t<<endl;
  //cout<<endl<<"updated ab matrix is:"<<endl;
  //display(ab1,row1,col2);
  //cout<<endl;

  for(int i=0;i<row;i++)
     for(int j=0;j<col;j++)
  {
      L[i][j]=x[i][j]-ab1[i][j];
  }

  //cout<<endl<<"updated L matrix is:"<<endl;
  //display(L,row,col);
  //cout<<endl;
 //cout<<endl<<"The dangerous value is:"<<endl;
 //double frobenius=frobeniusnorm(L,row,col);
// cout<<frobenius<<endl;


  }

}

/*3 3
4 3 1
1 2 3
4 5 6
3 1
6
7
8
1 3
3 6 7*/
