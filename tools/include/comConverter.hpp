#ifndef COM_CONVERTER_HPP
#define COM_CONVERTER_HPP

#include <vector>
#include <string>

class comConverter{
  public:
	comConverter();

	comConverter(const char *fname);

	bool load(const char *fname);

	double convertEject2lab(const double &com_);

	double convertRecoil2lab(const double &com_);

  private:
	std::vector<double> com;
	std::vector<double> ejectLab;
	std::vector<double> recoilLab;
	size_t length;
};

#endif
