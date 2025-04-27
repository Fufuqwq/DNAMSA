#include <bits/stdc++.h>
using namespace std;

const int INF = 1e9;
const int GAP = 2;
const int MISMATCH = 3;
const int MATCH = 0;

// 比较两个字符的得分（gap 用 '-' 表示）
int score(char a, char b)
{
    if (a == '-' || b == '-')
        return GAP;
    return (a == b) ? MATCH : MISMATCH;
}

int main()
{
    string A = "ACA";
    string B = "GAA";
    string C = "GGA";

    int n1 = A.size(), n2 = B.size(), n3 = C.size();

    // 三维DP数组
    vector<vector<vector<int>>> dp(n1 + 1, vector<vector<int>>(n2 + 1, vector<int>(n3 + 1, INF)));
    dp[0][0][0] = 0;
    vector<vector<vector<string>>> str(n1 + 1, vector<vector<string>>(n2 + 1, vector<string>(n3 + 1)));

    // 枚举所有状态
    for (int i = 0; i <= n1; ++i)
    {
        for (int j = 0; j <= n2; ++j)
        {
            for (int k = 0; k <= n3; ++k)
            {

                // 尝试所有合法 move（不全是 gap）
                for (int mask = 1; mask < 8; ++mask)
                {
                    int ni = i + ((mask & 1) ? 1 : 0);
                    int nj = j + ((mask & 2) ? 1 : 0);
                    int nk = k + ((mask & 4) ? 1 : 0);

                    if (ni > n1 || nj > n2 || nk > n3)
                        continue;

                    // 拿出当前字符或 gap
                    char ca = (mask & 1) ? A[i] : '-';
                    char cb = (mask & 2) ? B[j] : '-';
                    char cc = (mask & 4) ? C[k] : '-';

                    // 三两比对得分
                    int cost = score(ca, cb) + score(cb, cc) + score(ca, cc);

                    dp[ni][nj][nk] = min(dp[ni][nj][nk], dp[i][j][k] + cost);
                    if (dp[ni][nj][nk] == dp[i][j][k] + cost)
                    {
                        str[ni][nj][nk] = str[i][j][k] + (char)('0' + mask);
                    }
                }
            }
        }
    }

    cout << "最小比对代价: " << dp[n1][n2][n3] << '\n';
    cout << str[n1][n2][n3] << '\n';

    return 0;
}

/*
-ACA
GA-A
G-GA


*/