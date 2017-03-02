//#include <conio.h>
//#include <iostream.h>
#include <math.h>

//---------------------------------------------------
//	calculate minor of matrix OR build new matrix : k-had = minor

void minor2(float b[100][100],float a[100][100],int i,int n){
    int j,l,h=0,k=0;
    for(l=1;l<n;l++)
        for( j=0;j<n;j++){
            if(j == i)
                continue;
            b[h][k] = a[l][j];
            k++;
            if(k == (n-1)){
                h++;
                k=0;
            }
        }
}// end function

//---------------------------------------------------
//	calculate determinte of matrix
float det(float a[100][100],int n){
    int i;
    float b[100][100],sum=0;
    if (n == 1)
        return a[0][0];
    else if(n == 2)
        return (a[0][0]*a[1][1]-a[0][1]*a[1][0]);
    else
        for(i=0;i<n;i++){
            minor2(b,a,i,n);	// read function
            sum = (float) (sum+a[0][i]*pow(-1,i)*det(b,(n-1)));	// sum = determinte matrix
        }
    return sum;
}

//---------------------------------------------------
//	calculate transpose of matrix
void transpose(float c[100][100],float d[100][100],int n,float det){
    int i,j;
    float b[100][100];
    for (i=0;i<n;i++)
        for (j=0;j<n;j++)
            b[i][j] = c[j][i];
    for (i=0;i<n;i++)
        for (j=0;j<n;j++)
            d[i][j] = b[i][j]/det;	// array d[][] = inverse matrix
}// end function

//---------------------------------------------------
//	calculate cofactor of matrix
void cofactor(float a[100][100],float d[100][100],int n,float determinte){
    float b[100][100],c[100][100];
    int l,h,m,k,i,j;
    for (h=0;h<n;h++)
        for (l=0;l<n;l++){
            m=0;
            k=0;
            for (i=0;i<n;i++)
                for (j=0;j<n;j++)
                    if (i != h && j != l){
                        b[m][k]=a[i][j];
                        if (k<(n-2))
                            k++;
                        else{
                            k=0;
                            m++;
                        }
                    }
            c[h][l] = (float) pow(-1,(h+l))*det(b,(n-1));	// c = cofactor Matrix
        }
    transpose(c,d,n,determinte);	// read function
}// end function

//---------------------------------------------------
//	calculate inverse of matrix
std::vector<std::vector<double>> inverse(std::vector<std::vector<double>> &a)
{
    int n = a.size();

    float aa[100][100];
    float inv[100][100];

    for (int i = 0; i < n; i ++) {
        for (int j = 0; j < n; j ++) {
            aa[i][j] = a[i][j];
            inv[i][j] = 0;
        }
    }
    float deter = det(aa, n);
    if(deter == 0)
        inv[0][0] = 0;//   cout<<"\nInverse of Entered Matrix is not possible\n";
    else if(n == 1)
        inv[0][0] = 1;
    else
        cofactor(aa,inv,n,deter);	// read function

    std::vector<std::vector<double>> result(n, std::vector<double>(n, 0.0));
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            result[i][j] = inv[i][j];
        }
    }

    return result;
}

//---------------------------------------------------
//main fuction exe
/*
int main(void){
    int i,j,n;
    float a[100][100],d[100][100],deter;
    cout<<"\n C++ Program To Find Inverse Of Matrix\n\n";
    n = scanf(a);	// read function
    int print_matrix = 1;
    printf(a,n,print_matrix);	// read function
    deter = (float) det(a,n);	// read function
        cout<<"----------------------------------------------------\n";
        cout<<"\n\n Determinant of the Matrix : "<<deter;
        cout<<"\n\n-----------------------\n";
    inverse(a,d,n,deter);	// read function
    int print_inverse = 2;
    printf(d,n,print_inverse);	// read function
        cout<<"\n\n==============================* THE END *==============================\n";
        cout<<"\n		**** Thanks For Using The Program!!! ****\n";
    getch();
    return 0;
}// end main
*/
