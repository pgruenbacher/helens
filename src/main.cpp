#include <fstream>
#include <cxxopts.hpp>
#include <string>
#include <filesystem>
#include <iostream>
#include <cassert>
#include <cmath>

namespace fs = std::filesystem;

using InputPosType = int;

using HeightRaw = uint8_t;

class Terrain {
	std::vector<HeightRaw> m_data;
	int m_width{512};
	int m_height{512};
	float m_scale{11.f};
	float m_horizScale{30.f};
private:
	float getHeightFieldValue(uint32_t x, uint32_t y) const;
	HeightRaw getRawHeightFieldValue(uint32_t x, uint32_t y) const;
public:
	auto width() const { return m_width; }
	auto height() const { return m_height; }
	float horizScale() const { return m_horizScale; }
	void load(fs::path fpath);
	float interpolateHeight(float x, float y) const;
};



void Terrain::load(fs::path fpath) {
	std::ifstream fin(fpath, std::ios::binary);
	m_data.resize(m_width*m_height);
    // fin.open ("treeKc.raw", std::ios::in | std::ios::binary);
    while (!fin.eof())
    {
      fin.read( ( (char*)&m_data[0] ), m_data.size() * sizeof(m_data[0]) ) ;
      break;
    }
}

float lerp(float p_from, float p_to, float p_weight) { return p_from + (p_to - p_from) * p_weight; }

HeightRaw Terrain::getRawHeightFieldValue(uint32_t x, uint32_t y) const {
  assert(x < m_width);
  assert(y < m_height);
  return  m_data[(y * m_width) + x];
}

float Terrain::getHeightFieldValue(uint32_t x, uint32_t y) const {
	return getRawHeightFieldValue(x,y) * m_scale;
}

float Terrain::interpolateHeight(float x, float y) const {

// interpolate height bilinear
// https://www.pascalgamedevelopment.com/showthread.php?3512-interpolating-a-point-on-a-heightmap
// we'll probalby convert this to bicubic interpolation later....
// https://gamedev.stackexchange.com/questions/70322/how-do-i-avoid-interpolation-artefacts-when-scaling-up-a-heightmap


  x = std::clamp(x, 0.0f, m_width-1.0001f);
  y = std::clamp(y, 0.0f, m_height-1.0001f);

  // modf expects integral to be float but should be safe to cast to int after.
  float startX; // will be set by modf.
  float startY; // will be set by modf.
  float xFrac = modff(x, &startX);
  float yFrac = modff(y, &startY);


  assert(startX < m_width && startX >= 0);
  assert(startY < m_height && startY >= 0);
  uint32_t startXi = uint32_t(startX);
  uint32_t startYi = uint32_t(startY);

  // std::cout << "startXi " << startXi << " startYi " << startYi << " " << xFrac << " " << yFrac << std::endl;
  // test the border points for this
  assert(startXi <= (m_width-2));
  assert(startYi <= (m_height-2));

  auto h00 = getHeightFieldValue(startXi, startYi);
  auto h10 = getHeightFieldValue(startXi + 1, startYi);
  auto h01 = getHeightFieldValue(startXi, startYi + 1);
  auto h11 = getHeightFieldValue(startXi + 1, startYi + 1);

  return lerp(lerp(h00, h10, xFrac), lerp(h01, h11, xFrac), yFrac);
}

float sqr(float v) {
	return v*v;
}

float getHypotenuse(float xdiff, float ydiff) {
	// std::cout << "xdiff " << xdiff << " " << ydiff << std::endl;
	return sqrtf(sqr(xdiff) + sqr(ydiff));
}

float getDist(float x, float y, float x2, float y2) {
	// return sqrtf(std::pow(x-x2,2) + std::pow(y-y2,2));
	return getHypotenuse(x2-x, y2-y);
}

float getSurfaceDistance(float xa, float ya, float xb, float yb, const Terrain& terrain) {
	// std::cout << "xa " << xa << " " << xb << " " << ya << " " << yb << std::endl;
	const float dist = getDist(xa,ya,xb,yb);
	assert(dist > 0.f);
	const float numsteps = std::ceil(dist);
	assert(numsteps > 0);
	const float xderiv = (xb-xa)/numsteps;
	const float yderiv = (yb-ya)/numsteps;

	float totalSurfaceDist = 0.f;

	const float stepLength = sqrtf(std::pow(xderiv,2) + std::pow(yderiv,2)) * terrain.horizScale();
	// std::cout << "numsteps " << numsteps << std::endl;
	for (int step = 0; step < numsteps; ++step) {
		float x1 = xderiv*step + xa;
		float y1 = yderiv*step + ya;
		float x2 = xderiv*(step+1) + xa;
		float y2 = yderiv*(step+1) + ya;
		assert(x2 <= xb);
		assert(y2 <= yb);

		float h0 = terrain.interpolateHeight(x1,y1);
		float h1 = terrain.interpolateHeight(x2,y2);

		totalSurfaceDist += getHypotenuse(stepLength, h0-h1);
	}

	return totalSurfaceDist;
}

int main(int argc, const char** argv) {
	cxxopts::Options options("Helens", "parse and calculate surface distance between two points on two heightmaps");
	options.add_options()
		("pointa", "point A", cxxopts::value<std::vector<InputPosType>>())
		("pointb", "point B", cxxopts::value<std::vector<InputPosType>>())
		("premap", "heightmap raw data PRE", cxxopts::value<fs::path>()->default_value("../post.data"))
		("postmap", "heightmap raw data POST", cxxopts::value<fs::path>()->default_value("../pre.data"))
	;

	auto result = options.parse(argc, argv);


	if(!result.count("pointa")) {
		std::cerr << "--pointa required" << std::endl;
		return 1;
	}
	if(!result.count("pointb")) {
		std::cerr << "--pointb required" << std::endl;
		return 1;
	}

	auto pointa = result["pointa"].as<std::vector<InputPosType>>();
	auto pointb = result["pointb"].as<std::vector<InputPosType>>();

	if (pointa.size() != 2) {
		std::cerr << "pointa must be length 2" << std::endl;
		return 1;
	} 
	if (pointb.size() != 2) {
		std::cerr << "pointb must be length 2" << std::endl;
		return 1;
	} 

	std::cout << "pointA= " << pointa[0] << " " << pointa[1] << std::endl;
	std::cout << "pointB= " << pointb[0] << " " << pointb[1] << std::endl;

	auto premap = result["premap"].as<fs::path>();
	auto postmap = result["postmap"].as<fs::path>();

	assert(fs::exists(premap));
	assert(fs::exists(postmap));

	Terrain ta;
	ta.load(premap);

	if (pointa[0] < 0 || pointa[1] >= ta.width()) {
		std::cerr << "point a invalid " << std::endl;
		return 1;
	}
	if (pointb[0] < 0 || pointb[1] >= ta.height()) {
		std::cerr << "point b invalid " << std::endl;
		return 1;
	}

	float distA = getSurfaceDistance(pointa[0], pointa[1], pointb[0], pointb[1], ta);

	Terrain tb;
	tb.load(postmap);
	float distB = getSurfaceDistance(pointa[0], pointa[1], pointb[0], pointb[1], tb);

	std::cout << "Surface Dist of Terrain A = " << distA << std::endl;
	std::cout << "Surface Dist of Terrain B = " << distB << std::endl;
	std::cout << "Diff of terrains = " << (distA-distB) << std::endl;

	return 0;
}
