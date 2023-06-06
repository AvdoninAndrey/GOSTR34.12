#include "GOSTR3412.h"

uint8_t PoleGalua::findFirstOneBitInNumber(const int16_t& number, const uint8_t& startPositionFind)
{
	uint8_t find = 0;
	for (int j = startPositionFind; j >= 0; --j)
	{
		uint8_t dischargeTmp = (number >> j) & 0x1;
		if (dischargeTmp)
		{
			find = j;
			break;
		}
	}
	return find;
}

uint8_t PoleGalua::divisionWithModGaloiField(const uint16_t resultMultiplication)
{
	uint16_t resultMod = resultMultiplication;

	uint8_t currentStartPosition = 15;
	while (resultMod > 255)
	{
		for (int i = 7; i >= 0; --i)
		{
			bool flag = false;
			uint16_t tmp = multiplicationGaloisField<uint16_t>(this->IRREDUCIBLE_POLYNOMIAL, 1 << i);

			if (tmp >> currentStartPosition == 0)
			{
				uint8_t findFlag1 = findFirstOneBitInNumber(resultMod, currentStartPosition), findFlag2 = findFirstOneBitInNumber(tmp, currentStartPosition);

				if (findFlag1 == findFlag2)
				{
					resultMod ^= tmp;
					currentStartPosition = findFlag1;
					break;
				}
			}
		}
	}
	return static_cast<uint8_t>(resultMod);
}

void GOSTR3412::generationKey()
{
	/*
	первый полином - x^25 + x^3 + 1
	второй полином - x^127+ x^1 + 1
	*/
	srand(time(NULL));
	unsigned int firstRegister = 0; // степень полинома 25
	unsigned int secondRegister[4] = { 0, 0, 0, 0 }; // степень полинома 127

	firstRegister = rand() % 33554432; // начальное заполнение случайное число от 0 до 33554431 (33554431 в 2 ссч 25 единиц)
	if (firstRegister == 0)
	{
		firstRegister = 1;
	}

	for (int i = 3; i >= 0; --i) // начальное заполнение случайное число от 0 до 2 ^ 127 - 1  (будет храниться в массиве)
	{
		if (i == 0)
		{
			secondRegister[i] = rand() % 2147483648;
			if (secondRegister[i] == 0)
			{
				secondRegister[i] = 1;
			}
		}
		else
		{
			secondRegister[i] = rand() % 4294967296;
			if (secondRegister[i] == 0)
			{
				secondRegister[i] = 1;
			}
		}
	}

	for (int i = 0; i < this->SIZE_KEY / (sizeof(uint64_t) * 8); ++i)
	{
		uint64_t partOfKey = 0;
		for (int j = 0; j < sizeof(uint64_t) * 8; ++j)
		{
			unsigned int newBitFirstRegister = ((firstRegister >> 3) & 0x1) ^ (firstRegister & 0x1); // вычисляем бит обратной связи для первого регистра
			unsigned int newBitSecondRegister = ((secondRegister[3] >> 1) & 0x1) ^ (secondRegister[3] & 0x1); // вычисляем бит обратной связи для второго регистра
			
			unsigned int resultBit = newBitFirstRegister ^ newBitSecondRegister; // получаем результирующий бит

			partOfKey = (partOfKey << 1) | resultBit; // записываем полученный бит в часть ключ

			firstRegister = (firstRegister >> 1) | (resultBit << 24); // сдвигаем первый регистр

			for (int i = 0; i < 4; ++i)  // сдвигаем второй регистр
			{
				if (i == 0)
				{
					secondRegister[i] = (secondRegister[i] >> 1) | (resultBit << 30);
				}
				else
				{
					secondRegister[i] = (secondRegister[i] >> 1) | (secondRegister[i - 1] << 31);
				}
			}
		}
		this->key[i] = partOfKey;
	}
	std::ofstream keyFile("secret.key", std::ios::binary);
	for (int i = 0; i < 4; ++i)
	{
		keyFile.write(reinterpret_cast<char*>(&key[i]), 8);
	}
	keyFile.close();
}

uint64_t* GOSTR3412::conversionX(const uint64_t* inputVector, const uint64_t* gamma)
{
	uint64_t* outputVector = new uint64_t[2]{ 0 };
	for (int i = 0; i < 2; ++i)
	{
		outputVector[i] = inputVector[i] ^ gamma[i];
	}

	return outputVector;
}

uint64_t* GOSTR3412::conversionS(const uint64_t* inputVector)
{
	uint64_t* outputVector = new uint64_t[2]{ 0 };
	for (int i = 0; i < 2; ++i)
	{
		uint64_t partVector = inputVector[i];
		uint64_t partOutputVector = 0;
		for (int j = 7; j >= 0; --j)
		{
			partOutputVector = (partOutputVector << 8) | (kPi[(partVector >> (j * 8)) & 0xFF]);
		}
		outputVector[i] = partOutputVector;
	}
	return outputVector;
}

uint64_t* GOSTR3412::conversionR(const uint64_t* inputVector)
{
	uint64_t* outputVector = new uint64_t[2]{ 0 };
	PoleGalua operation;

	uint8_t summa = 0;
	int counter = 0;
	for (int i = 0; i < 2; ++i)
	{
		uint64_t partVector = inputVector[i];
		for (int k = 7; k >= 0; --k)
		{
			uint16_t tmp = operation.multiplicationGaloisField<uint8_t>(kL[counter], (partVector >> (k * 8) & 0xFF));
			if (tmp > 255)
			{
				summa ^= operation.divisionWithModGaloiField(tmp);
			}
			else
			{
				summa ^= static_cast<uint8_t>(tmp);
			}

			counter++;
		}
	}
	outputVector[0] = (static_cast<uint64_t>(summa) << 56) | (inputVector[0] >> 8);
	outputVector[1] = (inputVector[1] >> 8) | ((inputVector[0] & 0xFF) << 56);

	return outputVector;
}

uint64_t* GOSTR3412::conversionL(const uint64_t* inputVector)
{
	uint64_t* outputVector = new uint64_t[2]{ 0 };


	for (int i = 0; i < 16; ++i)
	{
		uint64_t* partVector;
		if (i == 0)
		{
			partVector = conversionR(inputVector);
		}
		else
		{
			partVector = conversionR(outputVector);
		}
		for (int j = 0; j < 2; ++j)
		{
			outputVector[j] = partVector[j];
		}
		delete[] partVector;
	}

	return outputVector;
}

GOSTR3412::GOSTR3412()
{
	generationConstant();
}

void GOSTR3412::generationConstant()
{
	for (int i = 1, indexIterationConstant = 0; i <= 32; ++i)
	{
		uint64_t firstPart = 0x0000000000000000, secondPart = i;

		uint64_t inputVectorI[2] = { firstPart , secondPart };
		uint64_t* resultConversionL = conversionL(inputVectorI);
		for (int j = 0; j < 2; ++j)
		{
			this->C[i - 1][j] = resultConversionL[j];
		}
		delete[] resultConversionL;
	}
}

void GOSTR3412::formationOfRoundKeys()
{
	uint64_t K1[2], K2[2] = { 0 };

	roundKeys[0][0] = K1[0] = key[0]; // первый раундовый ключ
	roundKeys[0][1] = K1[1] = key[1];

	roundKeys[1][0] = K2[0] = key[2]; // второй раундовый ключ
	roundKeys[1][1] = K2[1] = key[3];


	int indexRoundKey = 2;

	for (int i = 1; i <= 32; ++i)
	{
		uint64_t* resultL;
		uint64_t* resultXor;
		uint64_t* resultS;

		resultXor = conversionX(K1, C[i - 1]);
		resultS = conversionS(resultXor);
		resultL = conversionL(resultS);

		uint64_t* resultRoundFestel = conversionX(resultL, K2);


		for (int j = 0; j < 2; ++j)
		{
			K2[j] = K1[j];
			K1[j] = resultRoundFestel[j];
		}

		if (i % 8 == 0)
		{
			for (int k = 0; k < 2; ++k)
			{
				if (k == 0)
				{
					for (int j = 0; j < 2; ++j)
					{
						roundKeys[indexRoundKey][j] = K1[j];
					}
					indexRoundKey++;
				}
				else
				{
					for (int j = 0; j < 2; ++j)
					{
						roundKeys[indexRoundKey][j] = K2[j];
					}
					indexRoundKey++;
				}
			}
		}
		delete[] resultL, resultXor, resultS, resultRoundFestel;
	}
}


uint64_t* GOSTR3412::encryptBlockText(const uint64_t* blockPlainText)
{
	uint64_t* encryptResult = new uint64_t[2]{ 0 };
	uint64_t* resultRound = nullptr;
	for (int i = 1; i <= 10; ++i) // 10 раундов
	{
		uint64_t* resultXor;
		uint64_t* resultS;

		if (i == 10) // если последний раунд шифрования
		{
			resultXor = conversionX(resultRound, roundKeys[i - 1]);
			for (int j = 0; j < 2; ++j)
			{
				encryptResult[j] = resultXor[j];
			}

			delete[] resultXor;
			continue;
		}

		if (i == 1)
		{
			resultXor = conversionX(blockPlainText, roundKeys[i - 1]);
		}
		else
		{
			resultXor = conversionX(resultRound, roundKeys[i - 1]);
		}

		resultS = conversionS(resultXor);
		resultRound = conversionL(resultS);

		delete[] resultXor, resultS;
	}
	delete[] resultRound;

	return encryptResult;
}

uint64_t* GOSTR3412::conversionReverseS(const uint64_t* inputVector)
{
	uint64_t* outputVector = new uint64_t[2]{ 0 };
	for (int i = 0; i < 2; ++i)
	{
		uint64_t partVector = inputVector[i];
		uint64_t partOutputVector = 0;
		for (int j = 7; j >= 0; --j)
		{
			partOutputVector = (partOutputVector << 8) | (kPiReverse[(partVector >> (j * 8)) & 0xFF]);
		}
		outputVector[i] = partOutputVector;
	}
	return outputVector;
}

uint64_t* GOSTR3412::conversionReverseR(const uint64_t* inputVector)
{
	uint64_t* outputVector = new uint64_t[2]{ 0 };
	uint8_t summa = (inputVector[0] >> 56) & 0xFF;


	outputVector[0] = (inputVector[0] << 8) | (inputVector[1] >> 56);
	outputVector[1] = (inputVector[1] << 8) | summa;

	PoleGalua operation;

	int counter = 0;
	summa = 0;

	for (int i = 0; i < 2; ++i)
	{
		uint64_t partVector = outputVector[i];
		for (int k = 7; k >= 0; --k)
		{
			uint16_t tmp = operation.multiplicationGaloisField<uint8_t>(this->kL[counter], (partVector >> (k * 8) & 0xFF));
			if (tmp > 255)
			{
				summa ^= operation.divisionWithModGaloiField(tmp);
			}
			else
			{
				summa ^= static_cast<uint8_t>(tmp);
			}
			counter++;
		}
	}
	outputVector[1] = (outputVector[1] & (72057594037927935 << 8)) | summa;

	return outputVector;
}

uint64_t* GOSTR3412::conversionReverseL(const uint64_t* inputVector)
{
	uint64_t* outputVector = new uint64_t[2]{ 0 };

	for (int i = 0; i < 16; ++i)
	{
		uint64_t* partVector = nullptr;
		if (i == 0)
		{
			partVector = conversionReverseR(inputVector);
		}
		else
		{
			partVector = conversionReverseR(outputVector);

		}
		for (int j = 0; j < 2; ++j)
		{
			outputVector[j] = partVector[j];
		}
		delete[] partVector;
	}
	return outputVector;
}

uint64_t* GOSTR3412::decryptBlockText(const uint64_t* blockCipherText)
{
	uint64_t* decryptResult = new uint64_t[2]{ 0 };

	uint64_t* resultRound = nullptr;

	for (int i = 1, k = 10; i <= 10; ++i, --k) // 10 раундов
	{
		uint64_t* resultXor;
		uint64_t* resultReverseL;

		if (i == 10) // если последний раунд расшифрования
		{
			resultXor = conversionX(resultRound, roundKeys[k - 1]);
			for (int j = 0; j < 2; ++j)
			{
				decryptResult[j] = resultXor[j];
			}
			delete[] resultXor;
			continue;
		}

		if (i == 1)
		{
			resultXor = conversionX(blockCipherText, roundKeys[k - 1]);
		}
		else
		{
			resultXor = conversionX(resultRound, roundKeys[k - 1]);
		}

		resultReverseL = conversionReverseL(resultXor);
		resultRound = conversionReverseS(resultReverseL);

		delete[] resultXor, resultReverseL;
	}
	delete[] resultRound;

	return decryptResult;
}


std::string GOSTR3412::procedureOfAdditionTwo(const std::string& inputPlainText)
{
	std::string resultProcedure = inputPlainText;
	uint32_t countBytesInPlainText = inputPlainText.length();
	uint8_t howManyAddBytes = 16 - (countBytesInPlainText % 16);
	resultProcedure.push_back(static_cast<unsigned char>(128));
	if (howManyAddBytes == 0)
	{
		for (uint8_t i = 0; i < 15; i++)
		{
			resultProcedure.push_back(static_cast<unsigned char>(0));
		}
	}
	else
	{
		for (uint8_t i = 0; i < howManyAddBytes - 1; i++)
		{
			resultProcedure.push_back(static_cast<unsigned char>(0));
		}
	}
	return resultProcedure;
}

void GOSTR3412::encryptPlainTextModeECB(const std::string plainText, const std::string & filePathToSave)
{
	formationOfRoundKeys();
	std::string resultEncryptStr = "";
	std::string resultProcedure = procedureOfAdditionTwo(plainText);
	int countBlocks = resultProcedure.length() / 16;

	uint32_t indexByte = 0;
	for (int i = 0; i < countBlocks; ++i)
	{
		uint64_t* blockPlainText = new uint64_t[2]{ 0 };
		for (int j = 0; j < 16; ++j)
		{
			if (j < 8)
			{
				blockPlainText[0] = (blockPlainText[0] << 8) | static_cast<unsigned char>(resultProcedure[indexByte]);
				++indexByte;
			}
			else
			{
				blockPlainText[1] = (blockPlainText[1] << 8) | static_cast<unsigned char>(resultProcedure[indexByte]);
				++indexByte;
			}
		}
		uint64_t* resultEncrypt = encryptBlockText(blockPlainText);
		for (int j = 15, k = 7; j >= 0; --j, --k)
		{
			if (j == 7)
			{
				k = 7;
			}
			if (j >= 8)
			{
				resultEncryptStr += static_cast<char>((resultEncrypt[0] >> (k * 8)) & 0xFF);
			}
			else
			{
				resultEncryptStr += static_cast<char>((resultEncrypt[1] >> (k * 8)) & 0xFF);
			}
		}

		delete[] blockPlainText, resultEncrypt;
	}
	std::ofstream result(filePathToSave + "\\encrypt.enc", std::ios::binary | std::ios::out);
	
	result.write(resultEncryptStr.c_str(), resultEncryptStr.length());

	std::cout << "Your message has been successfully encrypted.\nThe result is presented in the file at the following path " << filePathToSave + "\\encrypt.enc" << std::endl;

}

void GOSTR3412::flipTheBytes(uint64_t & blockBytes)
{
	uint8_t* tmpBits = new uint8_t[8]{ 0 };
	for (int i = 7; i >= 0; --i)
	{
		tmpBits[i] = (blockBytes >> (i * 8)) & 0xFF;
	}
	blockBytes = 0;
	for (int i = 0; i < 8; ++i)
	{
		blockBytes = (blockBytes << 8) | tmpBits[i];
	}
}

std::string GOSTR3412::decryptCipherTextModeECB(const std::string & filePathCipherText, const std::string & filePathKey)
{

	std::string resultDecryptStr = "";

	std::ifstream keyFile(filePathKey, std::ios::binary);
	for (int i = 0; i < 4; ++i)
	{
		keyFile.read(reinterpret_cast<char*>(&key[i]), 8);
	}
	keyFile.close();
	formationOfRoundKeys();
	std::ifstream fileCipherText(filePathCipherText, std::ios::binary | std::ios::ate);

	int32_t fileSize = fileCipherText.tellg();
	fileCipherText.seekg(std::ios::beg);

	for (int i = 0; i < fileSize / 16; ++i)
	{
		uint64_t* blockCipherText = new uint64_t[2]{ 0 };

	
		fileCipherText.read(reinterpret_cast<char*>(&blockCipherText[0]), 8);
		fileCipherText.read(reinterpret_cast<char*>(&blockCipherText[1]), 8);

		flipTheBytes(blockCipherText[0]);
		flipTheBytes(blockCipherText[1]);


		uint64_t* resultDecryptBlock = decryptBlockText(blockCipherText);
		
		for (int j = 15, k = 7; j >= 0; --j, --k)
		{
			if (j == 7)
			{
				k = 7;
			}
			if (j >= 8)
			{
				resultDecryptStr += static_cast<char>((resultDecryptBlock[0] >> (k * 8)) & 0xFF);
			}
			else
			{
				resultDecryptStr += static_cast<char>((resultDecryptBlock[1] >> (k * 8)) & 0xFF);
			}
		}
		delete[] blockCipherText, resultDecryptBlock;
	}

	return procedureOfAdditionTwoReverse(resultDecryptStr);
}
std::string GOSTR3412::procedureOfAdditionTwoReverse(const std::string& resultDecryptText)
{
	std::string result = "";
	uint32_t countByteStr = resultDecryptText.length();
	uint32_t indexEnd = countByteStr - 1;
	for (int i = 0; i < 16; i++)
	{
		uint8_t tmp = (resultDecryptText[indexEnd] >> i * 8) & 0xFF;
		if (tmp > 0)
		{
			break;
		}
		indexEnd--;
	}

	for (int i = 0; i < indexEnd; ++i)
	{
		result += resultDecryptText[i];
	}

	return result;
}

void GOSTR3412::test()
{
	//encryptPlainTextModeECB("hello world");
	std::string res = decryptCipherTextModeECB("D:\\Programming\\C++\\Kuznechik\\encrypt.enc", "secret.key");
	std::cout << res << std::endl;
	std::ofstream a("res.txt", std::ios::binary);
	a.write(res.c_str(), res.length());
	std::cout << res << std::endl;

	
}

void GOSTR3412::testOperationAlgorithm()
{
	// Данный метод реализует проверку всех преобразований алгоритма. Тестовые данные были взяты с официальной документации

	/*Преобразование S
		S(ffeeddccbbaa99881122334455667700) = b66cd8887d38e8d77765aeea0c9a7efc,
		S(b66cd8887d38e8d77765aeea0c9a7efc) = 559d8dd7bd06cbfe7e7b262523280d39,
		S(559d8dd7bd06cbfe7e7b262523280d39) = 0c3322fed531e4630d80ef5c5a81c50b,
	*/
	uint64_t testConversionS[3][2] = { 0 };
	testConversionS[0][0] = 0xffeeddccbbaa9988;
	testConversionS[0][1] = 0x1122334455667700;
	testConversionS[1][0] = 0xb66cd8887d38e8d7;
	testConversionS[1][1] = 0x7765aeea0c9a7efc;
	testConversionS[2][0] = 0x559d8dd7bd06cbfe;
	testConversionS[2][1] = 0x7e7b262523280d39;

	std::cout << "Checking the S conversion:" << std::endl;
	for (int i = 0; i < 3; ++i)
	{
		std::cout << "Input vector: ";
		for (int j = 0; j < 2; ++j)
		{
			std::cout << std::hex << std::setfill('0') << std::setw(16) << testConversionS[i][j];
		}
		std::cout << "\t";
		uint64_t* resultS = conversionS(testConversionS[i]);
		std::cout << "Output vector: ";
		for (int j = 0; j < 2; ++j)
		{
			std::cout << std::hex << std::setfill('0') << std::setw(16) << resultS[j];
		}
		std::cout << std::endl;
		delete[] resultS;
	}
	std::cout << std::endl;
	/*Преобразование R
		R(00000000000000000000000000000100) = 94000000000000000000000000000001,
		R(94000000000000000000000000000001) = a5940000000000000000000000000000,
		R(а5940000000000000000000000000000) = 64а59400000000000000000000000000,
	*/
	uint64_t testConversionR[3][2] = { 0 };
	testConversionR[0][0] = 0x0000000000000000;
	testConversionR[0][1] = 0x0000000000000100;
	testConversionR[1][0] = 0x9400000000000000;
	testConversionR[1][1] = 0x0000000000000001;
	testConversionR[2][0] = 0xA594000000000000;
	testConversionR[2][1] = 0x0000000000000000;
	std::cout << "Checking the R conversion:" << std::endl;
	for (int i = 0; i < 3; ++i)
	{
		std::cout << "Input vector: ";
		for (int j = 0; j < 2; ++j)
		{
			std::cout << std::hex << std::setfill('0') << std::setw(16) << testConversionR[i][j];
		}
		std::cout << "\t";
		uint64_t* resultR = conversionR(testConversionR[i]);
		std::cout << "Output vector: ";
		for (int j = 0; j < 2; ++j)
		{
			std::cout << std::hex << std::setfill('0') << std::setw(16) << resultR[j];
		}
		std::cout << std::endl;
		delete[] resultR;
	}
	std::cout << std::endl;
	/*Преобразование L
		L(64a59400000000000000000000000000) = d456584dd0e3e84cc3166e4b7fa2890d,
		L(d456584dd0e3e84cc3166e4b7fa2890d) = 79d26221b87b584cd42fbc4ffea5de9a,
		L(79d26221b87b584cd42fbc4ffea5de9a) = 0e93691a0cfc60408b7b68f66b513c13,
	*/
	uint64_t testConversionL[3][2] = { 0 };
	testConversionL[0][0] = 0x64a5940000000000;
	testConversionL[0][1] = 0x0000000000000000;
	testConversionL[1][0] = 0xd456584dd0e3e84c;
	testConversionL[1][1] = 0xc3166e4b7fa2890d;
	testConversionL[2][0] = 0x79d26221b87b584c;
	testConversionL[2][1] = 0xd42fbc4ffea5de9a;

	std::cout << "Checking the L conversion:" << std::endl;
	for (int i = 0; i < 3; ++i)
	{
		std::cout << "Input vector: ";
		for (int j = 0; j < 2; ++j)
		{
			std::cout << std::hex << std::setfill('0') << std::setw(16) << testConversionL[i][j];
		}
		std::cout << "\t";
		uint64_t* resultL = conversionL(testConversionL[i]);
		std::cout << "Output vector: ";
		for (int j = 0; j < 2; ++j)
		{
			std::cout << std::hex << std::setfill('0') << std::setw(16) << resultL[j];
		}
		std::cout << std::endl;
		delete[] resultL;
	}
	std::cout << std::endl;
	/*Алгоритм развертывания ключа
	В качества тестового ключа использовался следующий ключ (256 бит): 8899aabbccddeeff 0011223344556677 fedcba9876543210 0123456789abcdef
		K1 = 8899aabbccddeeff0011223344556677,
		K2 = fedcba98765432100123456789abcdef,
		K3 = db31485315694343228d6aef8cc78c44,
		K4 = 3d4553d8e9cfec6815ebadc40a9ffd04,
		K5 = 57646468c44a5e28d3e59246f429f1ac,
		K6 = bd079435165c6432b532e82834da581b,
		K7 = 51e640757e8745de705727265a0098b1,
		K8 = 5a7925017b9fdd3ed72a91a22286f984,
		K9 = bb44e25378c73123a5f32f73cdb6e517,
		K10 = 72e9dd7416bcf45b755dbaa88e4a4043.
		Проверка правильности развёртывания ключей также проверяла выработку итерационных констант (С1-С32)
	*/
	this->key[0] = 0x8899aabbccddeeff;
	this->key[1] = 0x0011223344556677;
	this->key[2] = 0xfedcba9876543210;
	this->key[3] = 0x0123456789abcdef;
	generationConstant();
	formationOfRoundKeys();
	std::cout << "Checking the key deployment:" << std::endl;
	for (int i = 0; i < 10; ++i)
	{
		std::cout << "K" << i + 1 << " = ";
		for (int j = 0; j < 2; ++j)
		{
			std::cout << std::hex << std::setfill('0') << std::setw(16) << roundKeys[i][j];
		}
		std::cout << std::dec << std::endl;
	}
	std::cout << std::endl;
	/*Шифрование блока
	В качества тестового ключа использовался ключ из предыдущего пункта проверки.
	Открый текст = 1122334455667700 ffeeddccbbaa9988;
	Результат = 7f679d90bebc24305a468d42b9d4edcd.
	*/
	uint64_t blockPlainText[2] = { 0x1122334455667700, 0xffeeddccbbaa9988 };
	std::cout << "Checking encryption:" << std::endl;
	std::cout << "Plain text: ";
	for (int i = 0; i < 2; ++i)
	{
		std::cout << std::hex << std::setfill('0') << std::setw(16) << blockPlainText[i];
	}

	uint64_t* resultEncryption = encryptBlockText(blockPlainText);
	std::cout << "\t";
	std::cout << "Cipher text: ";
	for (int i = 0; i < 2; ++i)
	{
		std::cout << std::hex << std::setfill('0') << std::setw(16) << resultEncryption[i];
	}
	std::cout << std::endl;
	delete[] resultEncryption;


	/*Шифрование блока
	В качества тестового ключа использовался ключ из предыдущего пункта проверки.
	Шифртекст = 7f679d90bebc2430 5a468d42b9d4edcd.
	Открый текст = 1122334455667700 ffeeddccbbaa9988;
	*/
	uint64_t blockCipherText[2] = { 0x7f679d90bebc2430, 0x5a468d42b9d4edcd };
	std::cout << "Checking decryption:" << std::endl;
	std::cout << "Cipher text: ";
	for (int i = 0; i < 2; ++i)
	{
		std::cout << std::hex << std::setfill('0') << std::setw(16) << blockCipherText[i];
	}

	uint64_t* resultDecryption = decryptBlockText(blockCipherText);
	std::cout << "\t";
	std::cout << "Plain text: ";
	for (int i = 0; i < 2; ++i)
	{
		std::cout << std::hex << std::setfill('0') << std::setw(16) << resultDecryption[i];
	}
	std::cout << std::endl;
	delete[] resultDecryption;

}

