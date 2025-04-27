#include <bits/stdc++.h>
using namespace std;

const string ALPHABET = "ACGT-";
const double INF = 1e9 + 7;
const double GAP = 2;
const double PI = acos(-1);

// 封装一个打分矩阵
struct ScoreMatrix
{
    map<pair<char, char>, int> scores;
    int gapPenalty;

    ScoreMatrix(int gap = GAP) : gapPenalty(gap)
    {
        string alphabet = "ACGT";
        for (char a : alphabet)
            for (char b : alphabet)
                scores[{a, b}] = (a == b) ? 0 : 3;
    }

    int calcScore(char a, char b)
    {
        return (a == '-' || b == '-') ? gapPenalty : scores[{a, b}];
    }
};

struct Column
{
    map<char, double> freq;
};
const Column GAPCOLUMN = {{{'-', 1.0}}};

struct Profile
{
    vector<Column> columns;
    vector<string> finalSeqs;
    vector<int> order;
};

// 构建每次合并后产生的新的频率表
Profile buildProfile(const vector<string> seqs, const vector<int> ord)
{
    Profile prof;
    prof.finalSeqs = seqs;
    int rows = seqs.size(), cols = seqs[0].size();

    for (int i = 0; i < cols; i++)
    {
        Column col;
        map<char, int> buc;
        for (const auto &str : seqs)
            buc[str[i]]++;
        for (auto [c, v] : buc)
            col.freq[c] = 1.0 * v / rows;
        prof.columns.push_back(col);
    }
    prof.order = ord;

    return prof;
}

// 计算 Profile-Profile 得分
double calcScoreBetweenColumns(const Column a, const Column b, ScoreMatrix &rule)
{
    double score = 0.0;
    for (auto [ac, af] : a.freq)
        for (auto [bc, bf] : b.freq)
            score += af * bf * rule.calcScore(ac, bc);
    return score;
}

Profile mergeProfiles(Profile &lef, Profile &rig, ScoreMatrix &rule)
{
    int n = lef.columns.size(), m = rig.columns.size();
    vector<vector<double>> dp(n + 1, vector<double>(m + 1));
    vector<vector<int>> from(n + 1, vector<int>(m + 1));
    vector<int> ord = lef.order;
    ord.insert(ord.end(), rig.order.begin(), rig.order.end());

    // 初始化
    dp[0][0] = 0;
    for (int i = 1; i <= n; i++)
    {
        dp[i][0] = dp[i - 1][0] + calcScoreBetweenColumns(lef.columns[i - 1], GAPCOLUMN, rule);
        from[i][0] = 1;
    }
    for (int j = 1; j <= m; j++)
    {
        dp[0][j] = dp[0][j - 1] + calcScoreBetweenColumns(GAPCOLUMN, rig.columns[j - 1], rule);
        from[0][j] = 2;
    }

    for (int i = 1; i <= n; i++)
        for (int j = 1; j <= m; j++)
        {
            double match = dp[i - 1][j - 1] + calcScoreBetweenColumns(lef.columns[i - 1], rig.columns[j - 1], rule);
            double gapLef = dp[i][j - 1] + calcScoreBetweenColumns(GAPCOLUMN, rig.columns[j - 1], rule);
            double gapRig = dp[i - 1][j] + calcScoreBetweenColumns(lef.columns[i - 1], GAPCOLUMN, rule);

            if (match <= gapLef && match <= gapRig)
            {
                dp[i][j] = match;
                from[i][j] = 0;
            }
            else if (gapRig <= gapLef)
            {
                dp[i][j] = gapRig;
                from[i][j] = 1;
            }
            else
            {
                dp[i][j] = gapLef;
                from[i][j] = 2;
            }
        }

    int sizLef = lef.finalSeqs.size(), sizRig = rig.finalSeqs.size();
    vector<string> newSeqLef(sizLef, ""), newSeqRig(sizRig, "");

    int i = n, j = m;
    while (i > 0 || j > 0)
    {
        if (i > 0 && j > 0 && from[i][j] == 0)
        {
            --i, --j;
            for (int k = 0; k < sizLef; k++)
                newSeqLef[k] += lef.finalSeqs[k][i];
            for (int k = 0; k < sizRig; k++)
                newSeqRig[k] += rig.finalSeqs[k][j];
        }
        else if (i > 0 && from[i][j] == 1)
        {
            --i;
            for (int k = 0; k < sizLef; k++)
                newSeqLef[k] += lef.finalSeqs[k][i];
            for (int k = 0; k < sizRig; k++)
                newSeqRig[k] += '-';
        }
        else if (j > 0 && from[i][j] == 2)
        {
            --j;
            for (int k = 0; k < sizLef; k++)
                newSeqLef[k] += '-';
            for (int k = 0; k < sizRig; k++)
                newSeqRig[k] += rig.finalSeqs[k][j];
        }
        else
        {
            cerr << "Error In MergePorfiles\n";
        }
    }

    for (auto &str : newSeqLef)
        reverse(str.begin(), str.end());
    for (auto &str : newSeqRig)
        reverse(str.begin(), str.end());
    newSeqLef.insert(newSeqLef.end(), newSeqRig.begin(), newSeqRig.end());

    return buildProfile(newSeqLef, ord);
}

struct TreeNode
{
    Profile profile;
    int lson, rson;
    TreeNode() : lson(-1), rson(-1) {}
};

void postProgressiveAlign(int rt, vector<TreeNode> &tr, ScoreMatrix &rule)
{
    if (tr[rt].lson == -1 && tr[rt].rson == -1)
        return;
    int lef = tr[rt].lson, rig = tr[rt].rson;
    postProgressiveAlign(lef, tr, rule);
    postProgressiveAlign(rig, tr, rule);
    tr[rt].profile = mergeProfiles(tr[lef].profile, tr[rig].profile, rule);
}

void buildGuideTree(const vector<string> &seqs, vector<vector<double>> &distMatrix, vector<TreeNode> &tr)
{
    int n = seqs.size();
    int currentNum = n;
    vector<int> buc(n), siz(n + n - 1, 0);
    iota(buc.begin(), buc.end(), 0);
    fill(siz.begin(), siz.begin() + n, 1);

    for (int round = 1; round < n; ++round)
    {
        double minDist = INF;
        int lef = -1, rig = -1;
        for (auto id1 : buc)
            for (auto id2 : buc)
            {
                if (id1 == id2)
                    continue;
                if (distMatrix[id1][id2] < minDist)
                {
                    minDist = distMatrix[id1][id2];
                    lef = id1;
                    rig = id2;
                }
            }

        cout << lef << ' ' << rig << '\n';

        for (int id : buc)
        {
            if (id == lef || id == rig)
                continue;
            double disLef = distMatrix[id][lef];
            double disRig = distMatrix[id][rig];
            double disAvg = (disLef * siz[lef] + disRig * siz[rig]) / (siz[lef] + siz[rig]);
            distMatrix[id][currentNum] = distMatrix[currentNum][id] = disAvg;
        }

        tr[currentNum].lson = lef, tr[currentNum].rson = rig;
        siz[currentNum] = siz[lef] + siz[rig];

        buc.erase(find(buc.begin(), buc.end(), lef));
        buc.erase(find(buc.begin(), buc.end(), rig));

        buc.push_back(currentNum);

        ++currentNum;
    }

    assert(buc.size() == 1);
}

// Kmer+FTT加速模糊距离计算
void change(vector<complex<double>> &a, int n)
{
    vector<int> R(n);
    for (int i = 0; i < n; i++)
        R[i] = R[i / 2] / 2 + (i & 1) * n / 2;
    for (int i = 0; i < n; i++)
        if (i < R[i])
            swap(a[i], a[R[i]]);
}

void FFT(vector<complex<double>> &a, int n, int op)
{
    change(a, n);
    for (int m = 2; m <= n; m <<= 1)
    {
        complex<double> w1(cos(2 * PI / m), sin(2 * PI / m) * op);
        complex<double> x, y;
        for (int i = 0; i < n; i += m)
        {
            complex<double> wk(1, 0);
            for (int j = 0; j < m / 2; j++)
            {
                x = a[i + j];
                y = a[i + j + m / 2] * wk;
                a[i + j] = x + y;
                a[i + j + m / 2] = x - y;
                wk = wk * w1;
            }
        }
    }

    if (op == -1)
    {
        for (int i = 0; i < n; i++)
            a[i] /= n;
    }
}

int encodeChar(char c)
{
    switch (c)
    {
    case 'A':
        return 0;
    case 'C':
        return 1;
    case 'G':
        return 2;
    case 'T':
        return 3;
    default:
        cerr << "Error : Exist Char Out Of The ALPHABET\n";
        return 0;
    }
}

int kmerToIndex(const string &str, int st, int k)
{
    int id = 0;
    for (int i = 0; i < k; i++)
        id = id * 4 + encodeChar(str[st + i]);
    return id;
}

vector<complex<double>> encodeKmerFreq(const string &str, int k, int len)
{
    vector<complex<double>> vec(len, 0);
    for (int i = 0; i + k <= (int)str.size(); ++i)
    {
        int id = kmerToIndex(str, i, k);
        vec[id] += complex<double>(1, 0);
    }
    return vec;
}

double calcSimilaryDistFFT(const string &lef, const string &rig, int k = 3)
{
    int len = 1;
    while (len <= (1 << (2 * k)))
        len <<= 1;

    vector<complex<double>> vecLef = encodeKmerFreq(lef, k, len);
    vector<complex<double>> vecRig = encodeKmerFreq(rig, k, len);
    double dotProduct = 0.0, normL = 0.0, normR = 0.0;

    for (int i = 0; i < len; i++)
    {
        normL += norm(vecLef[i]);
        normR += norm(vecRig[i]);
    }

    FFT(vecLef, len, 1);
    FFT(vecRig, len, 1);
    for (int i = 0; i < len; i++)
        vecLef[i] *= vecRig[i];
    FFT(vecLef, len, -1);

    for (int i = 0; i < len; i++)
        dotProduct += vecLef[i].real();
    return 1.0 - dotProduct / (sqrt(normL) * sqrt(normR));
}

void initDistMatrix(const vector<string> &seqs, vector<vector<double>> &distMatrix)
{
    int n = seqs.size();
    for (int i = 0; i < n; i++)
        for (int j = i + 1; j < n; j++)
            distMatrix[i][j] = distMatrix[j][i] = calcSimilaryDistFFT(seqs[i], seqs[j], 1); // min(3,log(minLength)))
}

double calcDistBetweenNewSeq(const string &lef, const string &rig, ScoreMatrix &rule)
{
    double dist = 0.0;
    int n = lef.size();
    assert(lef.size() == rig.size());
    for (int i = 0; i < n; i++)
        dist += rule.calcScore(lef[i], rig[i]);
    return dist;
}

void reInitDistMatrix(const vector<string> &seqs, const vector<int> &order, vector<vector<double>> &distMatrix, ScoreMatrix &rule)
{
    int n = seqs.size();
    for (int i = 0; i < n; i++)
        for (int j = i + 1; j < n; j++)
        {
            int x = order[i], y = order[j];
            distMatrix[x][y] = distMatrix[y][x] = calcDistBetweenNewSeq(seqs[i], seqs[j], rule);
        }
}

double calculateScore(const vector<string> &seqs, ScoreMatrix &rule)
{
    double score = 0;
    int n = seqs[0].size(), k = seqs.size();
    for (int col = 0; col < n; ++col)
        for (int i = 0; i < k; i++)
            for (int j = i + 1; j < k; j++)
                score += rule.calcScore(seqs[i][col], seqs[j][col]);
    return score;
}

// 前序遍历进行拆边重构
void refine(int rt, Profile other, vector<TreeNode> &tr, ScoreMatrix &rule, Profile &finalProfile, double &finalScore)
{
    if (rt == -1)
        return;
    int lef = tr[rt].lson, rig = tr[rt].rson;
    if (lef == -1 && rig == -1) // 除了叶子节点外，别的节点都是有两个儿子的
        return;

    Profile newOther = mergeProfiles(tr[rig].profile, other, rule);
    Profile newAns = mergeProfiles(newOther, tr[lef].profile, rule);
    double newScore = calculateScore(newAns.finalSeqs, rule);
    if (newScore < finalScore)
    {
        finalProfile = newAns;
        finalScore = newScore;
    }
    refine(lef, newOther, tr, rule, finalProfile, finalScore);

    newOther = mergeProfiles(tr[lef].profile, other, rule);
    newAns = mergeProfiles(newOther, tr[rig].profile, rule);
    newScore = calculateScore(newAns.finalSeqs, rule);
    if (newScore < finalScore)
    {
        finalProfile = newAns;
        finalScore = newScore;
    }
    refine(rig, newOther, tr, rule, finalProfile, finalScore);
}

void solve()
{
    int n;
    cin >> n;
    vector<string> seqs(n, "");
    for (int i = 0; i < n; i++)
        cin >> seqs[i];
    vector<vector<double>> distMatrix(n + n - 1, vector<double>(n + n - 1, 0.0));
    vector<TreeNode> tr(n + n - 1);
    for (int i = 0; i < n; i++)
        tr[i].profile = buildProfile({seqs[i]}, {i});

    ScoreMatrix rule;
    initDistMatrix(seqs, distMatrix);
    buildGuideTree(seqs, distMatrix, tr);
    postProgressiveAlign(n + n - 2, tr, rule);

    Profile midProfile = tr.back().profile;
    vector<string> midSeqs = midProfile.finalSeqs;
    vector<int> midOrder = midProfile.order;

    for (int x : midOrder)
        cout << x << ' ';
    cout << '\n';

    reInitDistMatrix(midSeqs, midOrder, distMatrix, rule);
    buildGuideTree(seqs, distMatrix, tr);
    postProgressiveAlign(n + n - 2, tr, rule);

    Profile finalProfile = tr.back().profile;
    double finalScore = calculateScore(finalProfile.finalSeqs, rule);
    refine(n + n - 2, {}, tr, rule, finalProfile, finalScore);

    vector<string> finalSeqs = finalProfile.finalSeqs;

    for (string str : finalSeqs)
        cout << str << '\n';
    cout << finalScore << '\n';
}

#define fopen                                   \
    freopen("E:/vscode/oi/in.txt", "r", stdin); \
    freopen("E:/vscode/oi/out.txt", "w", stdout);
#define ios                  \
    ios::sync_with_stdio(0); \
    cin.tie(0);

signed main()
{

    fopen;

    solve();

    return 0;
}