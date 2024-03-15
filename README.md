# PACS-challenges
This repo contain the code used as solution of challenges given by the course PACS (Advanced Programming for Scientific Computing)

## How to use

I developed and test my code in the docker container provided by the course. To use the code, you need to have the docker installed in your machine. Then, you can pull the image using this command:
```bash 
docker pull albertoartoni1995/mk
```

Then, you can run the container using the following command:
```bash 
docker run --name pacs-env -v /path/to/host/folder:/home/jammy/shared-folder -it -d albertoartoni1995/mk
```

__TIPS & TRICKS:__
I have the issue that I can't compile or modify any file in the shared folder (I'm running the docker container in Arch Linux to have a well-set environment with all things done). To solve this, I run the following command to modify the permission of read/write of the folder:
```bash 
chmod -R 777 /path/to/host/folder
```
## Challenge 1

Develop the gradient descent method. 

More details in the README inside the folder challenge_1
