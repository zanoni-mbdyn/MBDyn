# MBDyn Tests Repository
This is a **private** repository for MBDyn test models. It is intended to contain all the
tests, either those that will have to remain private and those that can be made public. 

The procedure to make a test public is outlined in the next section. 

## How to publish a test (or a modification to an existing test)
Follow the steps below:
1. Make your modifications in the branch `master` of this repository; when you're done move to step (2.)
2. Checkout the `public` branch
   ```
   git checkout public
   ```
3. Checkout (do **not** merge!) your modifications from the `master` branch
   ```
   git checkout master <your_modifications>
   ```
4. Commit them
   ```
   git commit -m "Added public test"
   ```
5. Add the [mbdyn-tests-public](https://gitlab.com/zanoni-mbdyn/mbdyn-tests-public)
   repository as an additional remote (you can skip this step if `git remote -v` shows the
   public repository already)
   ```
   git remote add mbdyn-tests-public git@gitlab.com:zanoni-mbdyn/mbdyn-tests-public.git
   ```
6. Push the `public` branch of this repository to the `develop` branch of the
   `mbdyn-tests-public` repository
   ```
   git push mbdyn-tests-public public:develop
   ```
7. When you are done, in the `mbdyn-tests-public` repository, merge the `develop` branch into the `main`
   branch

