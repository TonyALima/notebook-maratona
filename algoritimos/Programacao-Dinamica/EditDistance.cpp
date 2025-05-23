#include <bits/stdc++.h>
using namespace std;

int min(int a, int b, int c)
{
    if (a < b && a < c)
        return a;
    if (b < c)
        return b;
    return c;
}

int editDistance(string s1, string s2)
{
    int m = s1.size();
    int n = s2.size();
    int PD[m + 1][n + 1];
    // Inicialização da primeira coluna (E(i, 0))
    for (int i = 0; i <= m; i++)
        PD[i][0] = i;
    // Inicialização da primeira linha (E(0, j))
    for (int j = 0; j <= n; j++)
        PD[0][j] = j;
    // Preenchendo a matriz dp
    for (int i = 1; i <= m; i++)
    {
        for (int j = 1; j <= n; j++)
        {
            if (s1[i - 1] == s2[j - 1])
            {
                PD[i][j] = PD[i - 1][j - 1]; // Sem custo extra
            }
            else
            {
                PD[i][j] = 1 + min(PD[i - 1][j], PD[i][j - 1], PD[i - 1][j - 1]);
            }
        }
    }
    return PD[m][n]; // Resposta final
}