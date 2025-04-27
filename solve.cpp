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

double calcDistBetweenSeq(const string &lef, const string &rig, ScoreMatrix &rule)
{
    int n = lef.size(), m = rig.size();
    vector<vector<double>> dp(n + 1, vector<double>(m + 1, 0.0));
    for (int i = 0; i <= n; i++)
        dp[i][0] = i * GAP;
    for (int j = 0; j <= m; j++)
        dp[0][j] = j * GAP;
    for (int i = 1; i <= n; i++)
        for (int j = 1; j <= m; j++)
        {
            double match = dp[i - 1][j - 1] + rule.calcScore(lef[i - 1], rig[j - 1]);
            double gapLef = dp[i][j - 1] + rule.calcScore('-', rig[j - 1]);
            double gapRig = dp[i - 1][j] + rule.calcScore(lef[i - 1], '-');
            dp[i][j] = min({match, gapLef, gapRig});
        }

    return dp[n][m];
}

void initDistMatrix(const vector<string> &seqs, vector<vector<double>> &distMatrix, ScoreMatrix &rule)
{
    int n = seqs.size();
    for (int i = 0; i < n; i++)
        for (int j = i + 1; j < n; j++)
            distMatrix[i][j] = distMatrix[j][i] = calcDistBetweenSeq(seqs[i], seqs[j], rule);
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

void refineByLeaveOneOut(const vector<string> finalSeqs, const vector<int> finalOrder, vector<string> &bestSeqs, double &bestScore, const vector<string> &seqs, ScoreMatrix &rule)
{
    int n = finalSeqs.size();
    for (int i = 0; i < n; i++)
    {
        string outOne = seqs[finalOrder[i]];
        vector<string> left;
        for (int j = 0; j < n; j++)
            if (j != i)
                left.push_back(finalSeqs[j]);
        Profile lef = buildProfile(left, {});
        Profile rig = buildProfile({outOne}, {});
        Profile newProfile = mergeProfiles(lef, rig, rule);
        double newScore = calculateScore(newProfile.finalSeqs, rule);
        if (newScore < bestScore)
        {
            bestSeqs = newProfile.finalSeqs;
            bestScore = newScore;
        }
    }
}

int getScore(char c, map<char, int> buc, ScoreMatrix &rule)
{
    int score = 0;
    for (auto [ch, num] : buc)
        score += (num - (ch == c ? 1 : 0)) * rule.calcScore(c, ch);
    return score;
}

void refineBySwapWithNeighborGap(vector<string> &finalSeqs, ScoreMatrix &rule)
{
    int n = finalSeqs.size(), m = finalSeqs[0].size();
    vector<vector<int>> ord(n);
    vector<map<char, int>> buc(m);

    for (int i = 0; i < n; i++)
    {
        ord[i].push_back(-1);
        for (int j = 0; j < m; j++)
            if (finalSeqs[i][j] != '-')
                ord[i].push_back(j);
        ord[i].push_back(m);
    }

    for (int col = 0; col < m; ++col)
        for (int row = 0; row < n; ++row)
            buc[col][finalSeqs[row][col]]++;

    while (1)
    {
        bool stable = true;
        for (int i = 0; i < n; i++)
            for (int j = 0; j < ord.size(); j++)
            {
                int p = ord[i][j];
                if (p == -1 || p == m)
                    continue;
                int l = ord[i][j - 1], r = ord[i][j + 1];
                for (int q = l + 1; q < r; ++q)
                {
                    int preScore = getScore(finalSeqs[i][p], buc[p], rule) + getScore(finalSeqs[i][q], buc[q], rule);
                    buc[p][finalSeqs[i][p]]--, buc[q][finalSeqs[i][q]]--;
                    swap(finalSeqs[i][q], finalSeqs[i][p]);
                    buc[p][finalSeqs[i][p]]++, buc[q][finalSeqs[i][q]]++;
                    int newScore = getScore(finalSeqs[i][p], buc[p], rule) + getScore(finalSeqs[i][q], buc[q], rule);
                    if (newScore < preScore)
                    {
                        stable = false;
                        ord[i][j] = q;
                        break;
                    }
                    else
                    {
                        buc[p][finalSeqs[i][p]]--, buc[q][finalSeqs[i][q]]--;
                        swap(finalSeqs[i][q], finalSeqs[i][p]);
                        buc[p][finalSeqs[i][p]]++, buc[q][finalSeqs[i][q]]++;
                    }
                }
            }

        if (stable)
            break;
    }
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
    initDistMatrix(seqs, distMatrix, rule);
    buildGuideTree(seqs, distMatrix, tr);
    postProgressiveAlign(n + n - 2, tr, rule);

    Profile finalProfile = tr.back().profile;
    vector<string> finalSeqs = finalProfile.finalSeqs;
    vector<int> finalOrder = finalProfile.order;

    for (string str : finalSeqs)
        cout << str << '\n';
    double finalScore = calculateScore(finalSeqs, rule);
    cout << finalScore << '\n';

    // refineByLeaveOneOut(finalSeqs, finalOrder, finalSeqs, finalScore, seqs, rule);

    refineBySwapWithNeighborGap(finalSeqs, rule);

    for (string str : finalSeqs)
        cout << str << '\n';
    finalScore = calculateScore(finalSeqs, rule);
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