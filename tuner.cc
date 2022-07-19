#include <bits/stdc++.h>
int auto_tuner(int NNZ, int N, int L, int *D)
{
    double x = (double)NNZ/N;
    double sum = 0;
    for (int l = 1; l <= L; l++)
    {
        sum = sum + D[l]*(l-x)*(l-x);
    }
    double sigma = sqrt(sum/N);
    int K = L;
    int *DD;
    DD = (int*)malloc(L * sizeof(int));
    DD[0]=0;
    for(int i = 2; i <= L; i++)
    {
        for(int j = 1; j<=L; j++)
        {
            DD[j]=D[j];
        }
        int count = 0;
        for(int j = i+1; j <= L; j++)
        {
            DD[i] = DD[i] + j/i*D[j];
            if(j%i == 0)
                count = count + (j/i-1)*D[j];
            else
            {
                count = count + (j/i)*D[j];
                DD[j%i] = DD[j%i] + D[j];
            }
        }
        double xx = NNZ/(N+count);
        double summ = 0;
        for(int l = 1; l <= i; l++)
        {
            summ = summ + DD[l]*(l-xx)*(l-xx);
        }
        double sigma2 = sqrt(summ/(N+count));
        if(sigma2 <= sigma/* && enough storage available*/) {
			sigma = sigma2;
			K = i;
		}
    }
	return K;
}


int main() {
	int D[] {0, 3, 1, 1};
	std::cout << auto_tuner(8, 3, 3, D) << std::endl;
	return 0;
}
