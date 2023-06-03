#include <iostream>
#include <filesystem>
#include <string>

bool createFolder(const std::string& folderPath) {
    // Check if the folder already exists
    if (std::filesystem::exists(folderPath)) {
        std::cout << "Folder already exists!" << std::endl;
        return false;
    }

    // Create the folder
    if (std::filesystem::create_directory(folderPath)) {
        std::cout << "Folder created successfully!" << std::endl;
        return true;
    }
    else {
        std::cout << "Failed to create folder!" << std::endl;
        return false;
    }
}

void copyHeaderFiles(const std::string& sourceDir, const std::string& destinationDir) {
    // Iterate over the files in the source directory
    for (const auto& entry : std::filesystem::directory_iterator(sourceDir)) {
        // Get the file path and name
        const std::string& filePath = entry.path().string();
        const std::string& fileName = entry.path().filename().string();

        // Check if the file ends with ".h" and starts with "PC"
        if (std::filesystem::is_regular_file(entry) && fileName.size() >= 2 && fileName.substr(0, 2) == "PC" && filePath.substr(filePath.size() - 2) == ".h") {
            // Generate the destination file path
            std::string destFilePath = destinationDir + "/" + fileName;

            // Copy the file to the destination directory
            std::filesystem::copy_file(filePath, destFilePath, std::filesystem::copy_options::overwrite_existing);

            std::cout << "Copied file: " << filePath << std::endl;
        }
    }
}

void copyDllFiles(const std::string& sourceDir, const std::string& destinationDir) {
    // Iterate over the files in the source directory
    for (const auto& entry : std::filesystem::directory_iterator(sourceDir)) {
        // Get the file path and name
        const std::string& filePath = entry.path().string();
        const std::string& fileName = entry.path().filename().string();

        // Check if the file is a DLL
        if (std::filesystem::is_regular_file(entry) && fileName.size() >= 4 && fileName.substr(fileName.size() - 4) == ".dll") {
            // Generate the destination file path
            std::string destFilePath = destinationDir + "/" + fileName;

            // Copy the file to the destination directory
            std::filesystem::copy_file(filePath, destFilePath, std::filesystem::copy_options::overwrite_existing);

            std::cout << "Copied file: " << filePath << std::endl;
        }
    }
}

bool renameFile(const std::string& filePath, const std::string& newFileName) {
    // Check if the file exists
    if (!std::filesystem::exists(filePath)) {
        std::cout << "File does not exist!" << std::endl;
        return false;
    }

    // Get the directory path of the file
    const std::string& directoryPath = std::filesystem::path(filePath).parent_path().string();

    // Generate the new file path
    std::string newFilePath = directoryPath + "/" + newFileName;

    // Rename the file
    std::filesystem::rename(filePath, newFilePath);

    std::cout << "File renamed successfully!" << std::endl;

    return true;
}


int main() {
    createFolder("../Build");
    createFolder("../Build/Practical-Components");
    createFolder("../Build/Practical-Components/include");
    createFolder("../Build/Practical-Components/include/PC");
    createFolder("../Build/Practical-Components/bin");
    createFolder("../Build/Practical-Components/bin/x64");
    createFolder("../Build/Practical-Components/bin/x64/Debug");
    createFolder("../Build/Practical-Components/bin/x64/Release");

    copyHeaderFiles("../PracticalComponents", "../Build/Practical-Components/include/PC");
    copyDllFiles("../x64/Debug", "../Build/Practical-Components/bin/x64/Debug");
    renameFile("../Build/Practical-Components/bin/x64/Debug/PracticalComponents.dll", "PracticalComponentsd.dll");

    copyDllFiles("../x64/Release", "../Build/Practical-Components/bin/x64/Release");
}


/*
folder structure
- Practical-Components/
    - include/
        - PC/
            - Header1.hpp
            - Header2.hpp
            - ...
    - bin/
        - Debug/
            - PracticalComponentsd.dll
        - Release/
            - PracticalComponents.dll
*/
