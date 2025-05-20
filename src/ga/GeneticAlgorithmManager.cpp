#include "GeneticAlgorithmManager.h"

GeneticAlgorithm::GeneticAlgorithm(int popSize_, int maxGen_, int tourSize_, double mutRate_, double s_, double s0_, double W_, double H_,const std::vector<double>& L_, const std::vector<double>& DY_, const std::vector<std::vector<double>>& hij_, const std::vector<double>& hi0_,  const std::vector<std::vector<double>>& Pij_, const std::vector<std::vector<double>>& Qij_, const std::vector<double>& de_)
{
	this->popSize = popSize_;
	this->maxGen = maxGen_;
	this->tourSize = tourSize;
	this->mutRate = mutRate_;
	this->s = s_;
	this->s0 = s0_;
	this->W = W_;
	this->H = H_;
	this->L = L_;
	this->DY = DY_;
	this->hij = hij_;
	this->hi0 = hi0_;
	this->Pij = Pij_;
	this->Qij = Qij_;
	this->de = de_;
	n = static_cast<int>(L.size());
	rng.seed(std::chrono::high_resolution_clock::now().time_since_epoch().count());
	initializePopulation();
}

GeneticAlgorithm::~GeneticAlgorithm()
{
}

void GeneticAlgorithm::initializePopulation()
{
	population.resize(popSize, std::vector<int>(n, 0));
	std::vector<int> base(n);
	for (int i = 0; i < n; ++i) 
		base[i] = i;
	for (int i = 0; i < popSize; ++i) 
	{
		shuffle(base.begin(), base.end(), rng);
		population[i] = base;
	}
}

int GeneticAlgorithm::tournamentSelection(const std::vector<double>& costs)
{
	std::uniform_int_distribution<int> distIdx(0, popSize - 1);
	int best = distIdx(rng);
	double bestCost = costs[best];
	for (int k = 1; k < tourSize; ++k) 
	{
		int cand = distIdx(rng);
		if (costs[cand] < bestCost) 
		{
			bestCost = costs[cand];
			best = cand;
		}
	}
	return best;
}

void GeneticAlgorithm::crossoverOX(const std::vector<int>& p1, const std::vector<int>& p2, std::vector<int>& c1, std::vector<int>& c2)
{
	std::uniform_int_distribution<int> distAB(0, n - 1);
	int a = distAB(rng);
	int b = distAB(rng);
	if (a > b) 
		std::swap(a, b);
	// 先复制父本区段
	for (int i = a; i <= b; ++i) 
	{
		c1[i] = p1[i];
		c2[i] = p2[i];
	}
	// 填充剩余基因
	auto fillChild = [&](const std::vector<int>& parent, std::vector<int>& child) {
		std::vector<bool> used(n, false);
		for (int i = a; i <= b; ++i) 
		{
			used[child[i]] = true;
		}
		int pos = (b + 1) % n;
		for (int i = 0; i < n; ++i) 
		{
			int gene = parent[(b + 1 + i) % n];
			if (!used[gene]) 
			{
				child[pos] = gene;
				pos = (pos + 1) % n;
			}
		}
		};
	fillChild(p2, c1);
	fillChild(p1, c2);
}

void GeneticAlgorithm::mutate(std::vector<int>& chromo)
{
	std::uniform_int_distribution<int> distIdx(0, n - 1);
	int i = distIdx(rng), j = distIdx(rng);
	std::swap(chromo[i], chromo[j]);
}

void GeneticAlgorithm::layoutPositions(const std::vector<int>& order, std::vector<std::pair<double, double>>& pos, int& rowCount, double& lastRowMaxWidth)
{
	pos.assign(this->n, { 0.0, 0.0 });
	std::vector<std::vector<int>> deviceInRow;  // 每一行的设备索引
	deviceInRow.reserve(n);

	int row = 0;
	deviceInRow.push_back(std::vector<int>{order[0]});
	// 第一台设备放置
	pos[order[0]].first = hi0[order[0]] + L[order[0]] / 2.0;  // x 坐标
	pos[order[0]].second = s0;                                // y 坐标
	double x_pos = pos[order[0]].first;
	int prev = order[0];

	// 依次放剩余设备
	for (int k = 1; k < n; ++k) 
	{
		int idx = order[k];		
		double next_x = x_pos + L[prev] / 2.0 + hij[prev][idx] + L[idx] / 2.0;           // 假设放在同一行
		if (next_x + L[idx] / 2.0 + de[idx] > W) 
		{
			row++;                           // 换行
			deviceInRow.push_back(std::vector<int>{idx});
			x_pos = hi0[idx] + L[idx] / 2.0;
		}
		else 
		{
			deviceInRow[row].push_back(idx);        // 继续同一行
			x_pos = next_x;
		}
		pos[idx].first = x_pos;
		pos[idx].second = s0 + row * s;
		prev = idx;
	}
	rowCount = row + 1; // 总行数

	// 计算最后一行设备中最大的宽度
	lastRowMaxWidth = 0.0;
	for (int idx : deviceInRow[row]) 
	{
		lastRowMaxWidth = std::max(lastRowMaxWidth, DY[idx]);
	}
}

double GeneticAlgorithm::fitness(const std::vector<int>& order)
{
	// 根据顺序计算坐标
	std::vector<std::pair<double, double>> posTmp;
	int m; double lsw;
	layoutPositions(order, posTmp, m, lsw);

	// 计算两两之间加权距离和
	double logisticsCost = 0.0;
	for (int i = 0; i < n; ++i) 
	{
		for (int j = i + 1; j < n; ++j) 
		{
			double dx = fabs(posTmp[i].first - posTmp[j].first);
			double dy = (posTmp[i].second == posTmp[j].second) ? 0.0: fabs(posTmp[i].second - posTmp[j].second);
			double dij = dx + dy;
			logisticsCost += Qij[i][j] * Pij[i][j] * dij;
		}
	}
	//  Y 方向越界惩罚
	double penalty = 0.0;
	if (s0 + (m - 1) * s + lsw > H) 
	{
		penalty = 1e9;
	}
	return logisticsCost + penalty;
}

void GeneticAlgorithm::run()
{
	std::vector<double> bestCostHistory(maxGen, 0.0);
	for (int gen = 0; gen < maxGen; ++gen) 
	{
		// 计算当前种群的所有个体成本
		std::vector<double> costs(popSize, 0.0);
		for (int i = 0; i < popSize; ++i)
		{
			costs[i] = fitness(population[i]);
		}
		// 记录这一代的最优成本
		double minCost = *min_element(costs.begin(), costs.end());
		bestCostHistory[gen] = minCost;

		//  生成新种群
		std::vector<std::vector<int>> newPop(popSize, std::vector<int>(n, 0));
		int idx = 0;
		while (idx < popSize) 
		{
			// 锦标赛选择两个父本
			int p1_idx = tournamentSelection(costs);
			int p2_idx = tournamentSelection(costs);
			const std::vector<int>& p1 = population[p1_idx];
			const std::vector<int>& p2 = population[p2_idx];

			//  OX 交叉产生两个子代
			std::vector<int> child1(n, -1), child2(n, -1);
			crossoverOX(p1, p2, child1, child2);

			// 变异
			if (std::uniform_real_distribution<double>(0.0, 1.0)(rng) < mutRate) 
			{
				mutate(child1);
			}
			if (std::uniform_real_distribution<double>(0.0, 1.0)(rng) < mutRate) 
			{
				mutate(child2);
			}

			// 放入新种群
			newPop[idx++] = child1;
			if (idx < popSize) 
			{
				newPop[idx++] = child2;
			}
		}
		population.swap(newPop);
	}

	// 最后一次迭代后，选出最佳个体并输出
	std::vector<double> finalCosts(popSize, 0.0);
	for (int i = 0; i < popSize; ++i) 
	{
		finalCosts[i] = fitness(population[i]);
	}
	int bestIdx = static_cast<int>(min_element(finalCosts.begin(), finalCosts.end()) - finalCosts.begin());
	bestOrder = population[bestIdx];

	// 计算最佳序列的布局坐标
	layoutPositions(bestOrder, bestPos, bestRowCount, bestLastRowMaxWidth);

	// 输出结果到屏幕（或写文件，根据需求自行修改）
	std::cout << "=== 最优设备排列顺序 (从 0 开始索引) ===" << std::endl;
	for (int id : bestOrder) 
	{
		std::cout << id << " ";
	}
	std::cout << "\n\n=== 最优布局坐标 (x, y) ===" << std::endl;
	for (int i = 0; i < n; ++i) 
	{
		std::cout << "设备 " << bestOrder[i] << ": ("
			<< bestPos[bestOrder[i]].first << ", "
			<< bestPos[bestOrder[i]].second << ")\n";
	}
}
