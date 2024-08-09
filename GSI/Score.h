#pragma once

#include <iostream>
#include <vector>

struct Identification {
	std::vector<std::size_t> peptides;
	const std::size_t spectrum;

	Identification(const std::size_t spectrum);
	Identification(const std::vector<std::size_t> peptides, const std::size_t spectrum);

	friend std::ostream& operator<<(std::ostream& os, const Identification& identification);

	void Add_Peptide(const std::size_t peptide);
};


struct Score : public Identification {
	double score;

	Score(const std::vector<std::size_t> peptides, const std::size_t spectrum ,double score);

	const Identification* Get_Edge() const;

	friend std::ostream& operator<<(std::ostream& os, const Score& score);
};

