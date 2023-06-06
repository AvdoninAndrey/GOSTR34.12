#include "GOSTR3412.h"



int main()
{
	setlocale(LC_ALL, "Russian");
	std::string enter;
	std::cout << "The program implements the encryption algorithm \"Grasshopper\" in the mode \"ECB\"\n";
	std::cout << "To check the operation of all transformations of the algorithm, enter \"test\"\n";
	std::cout << "To encrypt the information, you need to generate a key. Enter \"generation\" to generate the key\n";
	std::cout << "To decrypt, enter \"decrypt\"\n";
	std::cout << "to terminate the program, enter \"exit\"\n";
	while (true)
	{
		getline(std::cin, enter);
		if (enter == "generation")
		{
			GOSTR3412 encrypt;
			encrypt.generationKey();
			std::cout << "Enter the encryption message: ";
			std::string plainText, filePathToSave;
			getline(std::cin, plainText);
			std::cout << "Enter the directory to save the encryption result: ";
			getline(std::cin, filePathToSave);
			filePathToSave = "D:\\Programming\\C++\\Kuznechik";
			encrypt.encryptPlainTextModeECB(plainText, filePathToSave);
			std::cout << "Select an operation: ";
		}
		else if (enter == "test")
		{
			GOSTR3412 test;
			test.testOperationAlgorithm();
			std::cout << "Select an operation: ";

		}
		else if (enter == "decrypt")
		{
			GOSTR3412 decrypt;
			std::cout << "Specify the path to the file where the ciphertext is located: ";
			std::string filePath;
			getline(std::cin, filePath);
			std::cout << "specify the path to the file where the key is stored: ";
			std::string fileKeyPath;
			getline(std::cin, fileKeyPath);
			std::cout << "Result decrypt: " << decrypt.decryptCipherTextModeECB(filePath, fileKeyPath) << std::endl;;
			std::cout << "Select an operation: ";
		}
		else if (enter == "exit")
		{
			std::cout << "The program has finished its work\n";
			break;
		}
		else
		{
			std::cout << "Incorrect input, try again" << std::endl;
		}
	}

	return 0;
}
