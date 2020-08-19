# ConTaInR
  
## Description
ConTaInR is a Python based tool for taxonomic annotation of 16S rRNA sequences. The tool uses input data from the preprocessing section and runs BLAST against a 16S rRNA reference database. After running BLAST the results are analyzed and the taxonomic distribution is visualized in a sunburst diagram.<br>
Documentation is available on [Gitbook](https://minor2019.gitbook.io/minor2019/)

## Installation 

### Prerequisites
* Docker
* Docker-compose
* git

#### Download the repository
1. `git clone https://zero.han.nl/gitlab/Rogier.Stegeman/containr.git`
2. Switch to the right branch
   1. For production you need `master` so you're ready to go.
   2. For development switch to the current sprint e.g. `git checkout sprint-[sprint nr]`
      1. Then to create a new feature `git checkout [feature-name]` 
3. Copy the `.env.example` file and rename it to .env: `cp .env.example .env`. Fill in or change the variables to your preferences.

#### Usage
1.  Start the application: `docker-compose up`
2.  Go to the IP address 127.0.0.1:5000 or where the application is running to use the application.
2.  Stop the application: `Ctrl-C` and `docker-compose down`

## License
This project is licensed under the MIT License.

## Contact Information
*  Rogier Stegeman - [Github](https://github.com/RnRoger) - RD.Stegeman@student.han.nl
*  Alex Janse - [Github](https://github.com/grimcode) - a.janse@student.han.nl
*  Jordan Dekker
*  Robin van der Vliet