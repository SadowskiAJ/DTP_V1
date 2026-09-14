// Copyright (c) 2026, Dr Lijithan Kathirkamanathan and Dr Adam Jan Sadowski
// of Imperial College London and Dr Marc Seidel of Siemens Gamesa Renewable
// Energy. Developed with OpenAI Codex support.

// Copyright under a BSD 3-Clause License, see
// https://github.com/SadowskiAJ/DTP_V1.git

// Last modified at 09.38 on 14/09/2026


// A simple C++ program to transform the .pts text-based file containing tower coordinates
// into a .bin file with same but in a much more efficient (for storage and reading) binary format

// Compilation instructions (run from the directory containing this file):
//
// Windows, using a 64-bit MinGW-w64 installation:
//   g++ -O2 -std=c++17 -static -static-libgcc -static-libstdc++ PTS2BIN.cpp -o PTS2BIN.exe
//
// The static-linking options make the Windows executable independent of the
// compiler's libgcc and libstdc++ runtime DLLs.
//
// Linux, using GCC:
//   g++ -O2 -std=c++17 PTS2BIN.cpp -o PTS2BIN.exe
//
// The .exe filename is intentionally retained on Linux because the accompanying
// MATLAB function convertPTS2BIN.m currently invokes that exact filename. The
// resulting file is nevertheless a native Linux executable.


#include <cstdint>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>

int main(int argc, char** argv)
{
	if (argc != 2)
	{
		std::cerr << "Usage: PTS2BIN pts_file_root_no_extension" << std::endl;
		std::cerr << "e.g.: PTS2BIN test  (for test.pts)" << std::endl << std::endl;
		return 1;
	}

	const std::string sptsFileName = argv[1];
	std::ifstream ptsFile(sptsFileName + ".pts");
	if (!ptsFile.is_open())
	{
		std::cerr << "Unable to open input file " << sptsFileName + ".pts" << std::endl;
		return 2;
	}

	std::ofstream BINFile(sptsFileName + ".bin", std::ios::out | std::ios::binary | std::ios::trunc);
	if (!BINFile.is_open())
	{
		std::cerr << "Unable to create output file " << sptsFileName + ".bin" << std::endl;
		return 3;
	}

	std::cout << "File " << sptsFileName + ".pts" << " opened." << std::endl;
	std::uint64_t unCounter = 0;
	BINFile.write(reinterpret_cast<const char*>(&unCounter), sizeof(unCounter));

	double x = 0.0, y = 0.0, z = 0.0;
	std::string sLine;
	std::uint64_t lineNumber = 0;
	bool firstRecord = true;
	bool hasPointCount = false;
	std::uint64_t expectedPointCount = 0;
	while (std::getline(ptsFile, sLine))
	{
		++lineNumber;
		// Some exporters prefix the file with a UTF-8 byte-order mark.
		if (lineNumber == 1 && sLine.compare(0, 3, "\xEF\xBB\xBF") == 0)
		{
			sLine.erase(0, 3);
		}
		const auto first = sLine.find_first_not_of(" \t\r\n\v\f");
		if (first == std::string::npos)
		{
			continue;
		}

		// PTS may start with a single unsigned decimal point count. Only
		// recognise this on the first nonblank line; malformed rows must fail.
		if (firstRecord)
		{
			firstRecord = false;
			const auto last = sLine.find_last_not_of(" \t\r\n\v\f");
			const std::string record = sLine.substr(first, last - first + 1);
			if (record.find_first_not_of("0123456789") == std::string::npos)
			{
				hasPointCount = true;
				for (const char digit : record)
				{
					const auto value = static_cast<std::uint64_t>(digit - '0');
					if (expectedPointCount > (std::numeric_limits<std::uint64_t>::max() - value) / 10)
					{
						std::cerr << "Point count out of range on line " << lineNumber << std::endl;
						BINFile.close();
						std::remove((sptsFileName + ".bin").c_str());
						return 4;
					}
					expectedPointCount = expectedPointCount * 10 + value;
				}
				continue;
			}
		}

		std::istringstream isLine(sLine);
		if (!(isLine >> x >> y >> z))
		{
			std::cerr << "Invalid coordinate data on line " << lineNumber << std::endl;
			BINFile.close();
			std::remove((sptsFileName + ".bin").c_str());
			return 4;
		}

		BINFile.write(reinterpret_cast<const char*>(&x), sizeof(x));
		BINFile.write(reinterpret_cast<const char*>(&y), sizeof(y));
		BINFile.write(reinterpret_cast<const char*>(&z), sizeof(z));
		if (!BINFile)
		{
			std::cerr << "Write failure while creating " << sptsFileName + ".bin" << std::endl;
			BINFile.close();
			std::remove((sptsFileName + ".bin").c_str());
			return 5;
		}

		++unCounter;
		if (!(unCounter % 1000000))
		{
			std::cout << "Read " << unCounter / 1000000 << " million pts rows." << std::endl;
		}
	}

	if (ptsFile.bad() || (ptsFile.fail() && !ptsFile.eof()))
	{
		std::cerr << "Read failure in " << sptsFileName + ".pts" << std::endl;
		BINFile.close();
		std::remove((sptsFileName + ".bin").c_str());
		return 4;
	}
	if (hasPointCount && unCounter != expectedPointCount)
	{
		std::cerr << "Point count mismatch: header declares " << expectedPointCount
			<< " but read " << unCounter << " data points." << std::endl;
		BINFile.close();
		std::remove((sptsFileName + ".bin").c_str());
		return 4;
	}

	BINFile.seekp(0);
	BINFile.write(reinterpret_cast<const char*>(&unCounter), sizeof(unCounter));
	BINFile.close();
	if (!BINFile)
	{
		std::cerr << "Unable to finalise " << sptsFileName + ".bin" << std::endl;
		std::remove((sptsFileName + ".bin").c_str());
		return 6;
	}

	std::cout << "File " << sptsFileName + ".bin" << " created." << std::endl;
	std::cout << "Total " << unCounter << " data points." << std::endl << std::endl;

	return 0;
}
