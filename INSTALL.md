## Installation Requirements (Ubuntu and Mac):

* [R software environment](https://cran.r-project.org/)
* R packages
  * [shiny](https://shiny.rstudio.com/)
  * [shinyLP](https://github.com/jasdumas/shinyLP)
  * [phylocanvas](http://phylocanvas.org/)
  * [Kaphi](https://github.com/PoonLab/Kaphi)


## Requirements Installation Procedure (Ubuntu):

* The commands for each step are to be written/copied line by line to the terminal.

1. Updating and Upgrading The System  
    ```
    sudo apt-get update
    sudo apt-get upgrade
    ```
2. Installing R
    ```
    sudo apt-get install r-base
    sudo apt-get install r-base-dev
    ```
3. Installing R Packages
    ```
    R
    install.packages("shiny")
    install.packages("shinyLP")
    install.packages("phylocanvas")
    quit() 
    ```
4. Install [Kaphi](https://github.com/PoonLab/Kaphi) by following its [INSTALL.md](https://github.com/PoonLab/Kaphi/blob/master/INSTALL.md).


## Requirements Installation Procedure (Mac):

* The commands for each step are to be written/copied line by line to the terminal.

1. Install the latest version of R from the appropriate [mirror](https://cran.r-project.org/mirrors.html).
2. Installing Xcode
    ```
    xcode-select --install
    ```
   Follow the generated prompts to the end of installation. To verify that Xcode was correctly installed check what version    of Xcode was installed:
    ```
    xcodebuild -version
    ```
3. Installing Command Line Tools

   Go to http://developer.apple.com/downloads and sign in with your Apple ID (the same one you use for iTunes and app
   purchases). Search for "command line tools" (in the search field on the left), then click on version corresponding to the
   installed version of Xcode and click on the the .dmg link to download it. Run the .dmg and follow the generated prompts
   to the end of installation.
4. Installing Homebrew
    ```
    /usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"
    ```
5. Installing R Packages
    ```
    R
    install.packages("shiny")
    install.packages("shinyLP")
    install.packages("phylocanvas")
    quit() 
    ```
6. Install [Kaphi](https://github.com/PoonLab/Kaphi) by following its [INSTALL.md](https://github.com/PoonLab/Kaphi/blob/master/INSTALL.md).
    
    
## KaphiShiny Installation Procedure (Ubuntu and Mac):

* Navigate to your preferred location in the filesystem and clone KaphiShiny from the GitHhub repository
    ```
    git clone https://github.com/PoonLab/KaphiShiny
    ```
    
* Run KaphiShiny
    ```
    cd KaphiShiny
    Rscript -e "shiny::runApp(launch.browser=TRUE)"
    ```

