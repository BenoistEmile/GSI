#pragma once

#include <iostream>

struct Identification {
	const std::size_t peptide;
	const std::size_t spectrum;

	Identification(const std::size_t peptide, const std::size_t spectrum);

	friend std::ostream& operator<<(std::ostream& os, const Identification& identification);
};


struct Score : public Identification {
	double score;

	Score(const std::size_t peptide, const std::size_t spectrum ,double score);

	const Identification* Get_Edge() const;

	friend std::ostream& operator<<(std::ostream& os, const Score& score);
};

