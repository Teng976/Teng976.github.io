---
layout: post
title: Steps Creating a GitHub Pages site with Jekyll
---
1. install Ruby
2. gem install bundler
3. cd to-your-folder  
   bundle init
4. bundle config set --local path 'vendor/bundle'
5. bundle add jekyll
6. bundle exec jekyll new --force --skip-bundle .  
   bundle install 
7. update Gemfile  
   gem "github-pages", "~> VERSION", group: :jekyll_plugins
8. bundle exec jekyll serve  
   bundle exec jekyll serve --drafts // for displaying drafts
