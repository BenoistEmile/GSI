#pragma once

class JuliaSolver
{
  public:
  JuliaSolver();
  ~JuliaSolver();

  void addQVars(std::size_t size);
  void addXVars(std::size_t size);
  void addXConstraints(
    std::vector<std::vector<std::size_t> *> const &spectra_peptides);
  void addDeltaConstraints(
    std::vector<std::vector<std::tuple<std::size_t, std::size_t>> *> const
      &peptides_proteins,
    std::vector<float> const &useful_detectabilities,
    std::vector<std::vector<std::size_t> *> const &peptides_spectra,
    std::vector<unsigned int> const &peptides_sure);
  void addObjectiveDelta(
    std::vector<std::vector<std::tuple<std::size_t, std::size_t>> *> const
      &peptides_proteins,
    std::vector<float> const &useful_detectabilities,
    std::vector<std::vector<std::size_t> *> const &peptides_spectra,
    std::vector<unsigned int> const &peptides_sure,
    std::vector<std::vector<std::tuple<std::size_t, std::size_t>> *>
      const &useless_peptides_proteins,
    std::vector<float>  const &useless_detectabilities);
  void addObjectivePSMs(std::vector<const Score *> const &useful_scores);
  void addObjectiveParcimony();

    private: 
  std::size_t m_num_Q;
  std::size_t m_num_X;
};