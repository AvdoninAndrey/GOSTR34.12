Реализация блочного шифра ГОСТ Р 3212 "Кузнечик" в режиме "ECB";
====Зашифрование====
При зашифровании информации происходит генерация ключа (256 бит) на основе комбинирующего генератора ПСП.
Данный ключ записывается в файл с названием "secret.key".
После этого необходимо ввести сообщение (открытый текст) и каталог, в котором будет сохранён результат зашифрования.
Файл результата зашифрования имеет название "encrypt.enc".
====Расшифрование====
При расшифровании необходимо указать путь до файла с ключом и путь до файла с шифртекстом.
Результатом расшифрования является строка открытого текста.