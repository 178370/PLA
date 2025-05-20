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
	// �ȸ��Ƹ�������
	for (int i = a; i <= b; ++i) 
	{
		c1[i] = p1[i];
		c2[i] = p2[i];
	}
	// ���ʣ�����
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
	std::vector<std::vector<int>> deviceInRow;  // ÿһ�е��豸����
	deviceInRow.reserve(n);

	int row = 0;
	deviceInRow.push_back(std::vector<int>{order[0]});
	// ��һ̨�豸����
	pos[order[0]].first = hi0[order[0]] + L[order[0]] / 2.0;  // x ����
	pos[order[0]].second = s0;                                // y ����
	double x_pos = pos[order[0]].first;
	int prev = order[0];

	// ���η�ʣ���豸
	for (int k = 1; k < n; ++k) 
	{
		int idx = order[k];		
		double next_x = x_pos + L[prev] / 2.0 + hij[prev][idx] + L[idx] / 2.0;           // �������ͬһ��
		if (next_x + L[idx] / 2.0 + de[idx] > W) 
		{
			row++;                           // ����
			deviceInRow.push_back(std::vector<int>{idx});
			x_pos = hi0[idx] + L[idx] / 2.0;
		}
		else 
		{
			deviceInRow[row].push_back(idx);        // ����ͬһ��
			x_pos = next_x;
		}
		pos[idx].first = x_pos;
		pos[idx].second = s0 + row * s;
		prev = idx;
	}
	rowCount = row + 1; // ������

	// �������һ���豸�����Ŀ��
	lastRowMaxWidth = 0.0;
	for (int idx : deviceInRow[row]) 
	{
		lastRowMaxWidth = std::max(lastRowMaxWidth, DY[idx]);
	}
}

double GeneticAlgorithm::fitness(const std::vector<int>& order)
{
	// ����˳���������
	std::vector<std::pair<double, double>> posTmp;
	int m; double lsw;
	layoutPositions(order, posTmp, m, lsw);

	// ��������֮���Ȩ�����
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
	//  Y ����Խ��ͷ�
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
		// ���㵱ǰ��Ⱥ�����и���ɱ�
		std::vector<double> costs(popSize, 0.0);
		for (int i = 0; i < popSize; ++i)
		{
			costs[i] = fitness(population[i]);
		}
		// ��¼��һ�������ųɱ�
		double minCost = *min_element(costs.begin(), costs.end());
		bestCostHistory[gen] = minCost;

		//  ��������Ⱥ
		std::vector<std::vector<int>> newPop(popSize, std::vector<int>(n, 0));
		int idx = 0;
		while (idx < popSize) 
		{
			// ������ѡ����������
			int p1_idx = tournamentSelection(costs);
			int p2_idx = tournamentSelection(costs);
			const std::vector<int>& p1 = population[p1_idx];
			const std::vector<int>& p2 = population[p2_idx];

			//  OX ������������Ӵ�
			std::vector<int> child1(n, -1), child2(n, -1);
			crossoverOX(p1, p2, child1, child2);

			// ����
			if (std::uniform_real_distribution<double>(0.0, 1.0)(rng) < mutRate) 
			{
				mutate(child1);
			}
			if (std::uniform_real_distribution<double>(0.0, 1.0)(rng) < mutRate) 
			{
				mutate(child2);
			}

			// ��������Ⱥ
			newPop[idx++] = child1;
			if (idx < popSize) 
			{
				newPop[idx++] = child2;
			}
		}
		population.swap(newPop);
	}

	// ���һ�ε�����ѡ����Ѹ��岢���
	std::vector<double> finalCosts(popSize, 0.0);
	for (int i = 0; i < popSize; ++i) 
	{
		finalCosts[i] = fitness(population[i]);
	}
	int bestIdx = static_cast<int>(min_element(finalCosts.begin(), finalCosts.end()) - finalCosts.begin());
	bestOrder = population[bestIdx];

	// ����������еĲ�������
	layoutPositions(bestOrder, bestPos, bestRowCount, bestLastRowMaxWidth);

	// ����������Ļ����д�ļ����������������޸ģ�
	std::cout << "=== �����豸����˳�� (�� 0 ��ʼ����) ===" << std::endl;
	for (int id : bestOrder) 
	{
		std::cout << id << " ";
	}
	std::cout << "\n\n=== ���Ų������� (x, y) ===" << std::endl;
	for (int i = 0; i < n; ++i) 
	{
		std::cout << "�豸 " << bestOrder[i] << ": ("
			<< bestPos[bestOrder[i]].first << ", "
			<< bestPos[bestOrder[i]].second << ")\n";
	}
}
